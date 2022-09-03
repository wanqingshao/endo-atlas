import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.graph_objects as go
import dash_table as dt
import pandas as pd
import numpy as np
import json
import sqlite3
import textwrap


app = dash.Dash(__name__, external_stylesheets=[dbc.themes.CYBORG])
application = app.server



'''
~~~~~~~~~~~~~~~~~~~
~~ APP Data ~~
~~~~~~~~~~~~~~~~~~~
'''
def obj_data_to_mesh3d(odata):
    # odata is the string read from an obj file
    vertices = []
    faces = []
    #lines = odata.splitlines()

    for line in odata:
        slist = line.split()
        if slist:
            if slist[0] == 'v':
                vertex = np.array(slist[1:], dtype=float)
                vertices.append(vertex)
            elif slist[0] == 'f':
                face = []
                for k in range(1, len(slist)):
                    face.append([int(s) for s in slist[k].replace('//','/').split('/')])
                if len(face) > 3: # triangulate the n-polyonal face, n>3
                    faces.extend([[face[0][0]-1, face[k][0]-1, face[k+1][0]-1] for k in range(1, len(face)-1)])
                else:
                    faces.append([face[j][0]-1 for j in range(len(face))])
            else: pass


    return np.array(vertices), np.array(faces)

def get_verticies(odata):
    vertices = []
    for line in odata:
        slist = line.split()
        if slist:
            if slist[0] == 'v':
                vertex = np.array(slist[1:], dtype=float)
                vertices.append(vertex)
            else: pass
    return np.array(vertices)


obj_file = open("obj/segmentation grid_no_top_dot.obj", "r")
obj_data= obj_file.readlines()
vertices, faces = obj_data_to_mesh3d(obj_data)
x, y, z = vertices[:,:3].T
I, J, K = faces.T
x_pos = z
y_pos = x
z_pos = y

no_express = [0 for i in z_pos]

obj_file2 = open("obj/segmentation grid_grid_only3.obj", "r")
obj_data2= obj_file2.readlines()
vertices2 = get_verticies(obj_data2)
point_x, point_y, point_z = vertices2[:,:3].T


insitu_file = open("obj/insitu_list.txt", "r").readlines() ## need to be updated
insitu_gene = [x[:(len(x)-1)] for x in insitu_file]

'''
~~~~~~~~~~~~~~~~~~~
~~ APP Figures ~~
~~~~~~~~~~~~~~~~~~~
'''


defined_colorscale = "viridis"

with open('obj/20220116_vertex_dict.json', 'r') as f:
    vertex_dic = json.load(f)

with open('obj/20220130_id_to_name.json', 'r') as f:
    id_name_dic = json.load(f)

cor_df = pd.read_csv("obj/20210120_predicted_sig_top10_correlation_by_bin_51markers.csv")


prediction_db  = sqlite3.connect("obj/prediction_20220130_3d_51markers.db")
all_gene_names = id_name_dic.keys()
gene_id = "NV2g017667000.1"
db_list =  list(prediction_db.execute("SELECT GENENAME,PREDIC from PREDDICTION WHERE ID=?", (gene_id, )))

similar_gene_p = cor_df.loc[cor_df["Var1"] == gene_id, "top_genes"].values


def get_apple_plot(db_list, gene_id, defined_colorscale):

    mesh = go.Mesh3d(
        x=x_pos,y=y_pos,z=z_pos,
        vertexcolor=vertices[:, 3:], #the color codes must be triplets of floats  in [0,1]!!
        i=K,j=I,k=J,
        name='',
        showscale=True,
        lighting=dict(ambient= 0.18,
                      diffuse= 1,
                      fresnel=  .1,
                      specular= 1,
                      roughness= .4),
        lightposition=dict(x=100,y=200,z=150),
        colorscale=defined_colorscale,
        intensitymode='vertex',
        colorbar=dict(xanchor="left", x  = -0.1,
                      yanchor = "top", y = 1,
                      xpad = 0.5, len = 0.3, thickness = 4,
                      title = "Expression prediction")
        )
    layout = go.Layout(
        font=dict(size=14, color='white'),
        width=600,
        height=450,
        scene=dict(xaxis=dict(visible=False),
                   yaxis=dict(visible=False),
                   zaxis=dict(visible=False),
                   aspectratio=dict(x=1,y=1,z=0.8
                             ),
                   camera=dict(eye=dict(x=0, y=1, z=1.2))),
        margin=dict(t=50))

    try:
        gene_name = db_list[0][0]
        intensity = [float(x) for x in db_list[0][-1].split(",")]
        mesh.update(intensity = np.array(intensity))
    except:
        try:
            gene_name =id_name_dic[gene_id] + "<br>is not expressed in this study"
            mesh.update(intensity = no_express)
        except:
            gene_name = "Invalid gene name"
            mesh.update(intensity=no_express)

    newfig = go.Figure(data=[mesh], layout=layout)
    newfig.add_scatter3d(x=point_z, y=point_x,z = point_y,mode='markers',
                         marker=dict(size=1.5,color="black"))

    newfig.update_layout(
        title=gene_name,
        title_x = 0.5,
        paper_bgcolor='rgb(0,0,0)',
        font=dict(size=14, color='white'),
        hovermode = False,
        width = 400,
        margin=dict(l=40, r=40, t=50, b=20),
        height = 400)
    return newfig


def get_similar_genes(gene_id, cor_df, id_name_dic):
    similar_gene_df = cor_df.loc[cor_df["Var1"] == gene_id, :]
    similar_gene_df.columns = ["Query GeneID", "Similar GeneID", "Correlation coefficient"]
    similar_gene_df["Query GeneName"] = id_name_dic[gene_id]
    similar_gene_df["Similar GeneName"] = [id_name_dic[x] for x in similar_gene_df["Similar GeneID"]]
    reordered_name = ["Query GeneName",  "Similar GeneName",  "Correlation coefficient"]
    data = similar_gene_df.loc[:,reordered_name ].to_dict('records')
    return data

newfig = get_apple_plot(db_list, gene_id, defined_colorscale)

if( gene_id not in insitu_gene):
    insitu_url = 'assets/insitu_NA.png'
else:
    insitu_url = 'assets/' + gene_id + ".png"


similar_gene_df = get_similar_genes(gene_id, cor_df, id_name_dic)

'''
~~~~~~~~~~~~~~~~
~~ APP LAYOUT ~~
~~~~~~~~~~~~~~~~
'''

title = html.Div(
    [
    html.Div(
        children=[html.Img(src='https://raw.githubusercontent.com/wanqingshao/pics/master/logo.png', style={'height':'80pt', 'width':'80pt'})],
        style = {'float':'left'}),
    html.H2("Endo-Atlas"),
    html.H6("Spatial reconstruction of single-cell transcriptomic data enables systematic dissection of Hox dependent tissue segmentation in the starlet sea anemone Nematostella vectensis"),
    html.Br()
    ]
)

web_des = html.Div(
    html.Div(
        children=[
        html.H5("Project description"),
        dbc.Row(
           children = [
               html.Div(
                   html.Img(src='assets/intro.png', style={ 'width':'35vw'}),
                   style ={'margin': 5}),
           ],
           className="h-75",
           justify="center"
        ),
        html.P("The development of metameric body segments was central to the evolution of complex animal body plans. In bilaterian organisms, the capacity for individual segments to polarize allows for the further specialization of metameric structures in coordination with axial pattern. Precisely how segmentation and segment polarization pathways arose in evolution remains poorly understood. Here we investigate the ancient origins of segment polarization in the pre-bilaterian sea anemone Nematostella vectensis, where select endomesodermal structures display polarized pattern formation along the directive axis. To analyze endomesodermal gene expression on a global scale, we first used single cell RNA-seq data to construct a comprehensive 3-D gene expression atlas of the developing Nematostella endoderm. We then systematically validated the accuracy of this in silico prediction using fluorescent in situ hybridization and identified two homeobox genes previously implicated in bilaterian segmental polarization, Lbx and Uncx, that form reciprocal stripes demarcating each pair of larval segment boundaries.  Genetic perturbation of Lbx function led to the duplication of Uncx-side segment identity and resulted in a mirrored pattern of retractor muscles in developing polyps, indicating a key role for segmental polarity in development of the polarized pattern of muscles. Taken in sum, these experiments provide a global endomesodermal gene expression atlas and demonstrate the molecular basis for segment polarity program in a cnidarian organism, hinting that polarized segments were likely present in the Cnidaria-Bilateria common ancestor over 600 million years ago."),
        ],
        style={'marginRight' : 25}

    ),
    style={ 'border-color': "grey",'border-right-style': 'outset'})

dropdown = dcc.Dropdown(
            id = 'gene_list',
            options = [{'label':id_name_dic[key], 'value':key} for key in id_name_dic.keys()],
            style = {'width':'200pt'},
            value=gene_id
        )

apple_plot = html.Div(
    [

    html.Div(
        [
        html.H5("Endo-Atlas database"),
        ],
        style = {'textAlign':'center'}),

    html.Div(
        [
        dbc.Row(
            [
            dbc.Col(html.P("Please type or select a gene "), width={"offset": 3}),
            dbc.Col(dropdown),
            ],
            justify = "center"
        ),],
        style = {'marginTop': 15}
    ),
    dbc.Row([
        dbc.Col(
            [
            dbc.Row(html.H5("Spacial prediction"), style = {'marginLeft':15, "textAlign" : "center"}, justify = "center"),
            html.Br(),

            dbc.Row(
                dcc.Graph(
                    id='apple_plot',
                    figure= newfig,
                    config={'editable': True, 'scrollZoom': True},
                    style = {"align":'center'}),
                    justify = "center"
                )

            ],
        ),
        dbc.Col(
            [
            html.H5("In situ results"),
            html.Br(),
            html.Br(),
            html.H5(children=gene_id, id = "gene_name"),
            html.Img(
                id = "insitu",
                src=insitu_url,
                alt = "Expression pattern has not been experimentally validated",
                style={'height':'200pt', 'width':'200pt'})
            ],
            style = { 'marginLeft':15, 'textAlign':'center'}
        )
        ],
        justify="center"),

    html.Div(
        [
            html.H6("Top 10 genes with similar spatial pattern based on pattern correlation analysis"),
            html.Div(dt.DataTable(
                id='tbl', data=similar_gene_df,
                columns=[{'name': 'Query GeneName', 'id': 'Query GeneName'},
                         {'name': 'Similar GeneName', 'id': 'Similar GeneName'},
                         {'name': 'Correlation coefficient', 'id': 'Correlation coefficient'}],
                style_data={
                    'color': 'white',
                    'backgroundColor': 'black',
                    'whiteSpace': 'normal',
                    'height': 'auto',
                    "maxwidth": "20"
                },
                style_cell={
                    "maxwidth": "20",
                    'whiteSpace': 'normal',
                    'height': 'auto',
                    'textOverflow': 'ellipsis'
                },
                style_header={
                        'backgroundColor': 'rgb(71, 71, 71)',
                        'color': 'white',
                        'fontWeight': 'bold'
                    }
            ), style = {"marginRight":30}),

        ],
        style={'textAlign': 'center', 'marginTop':10 , "marginLeft":10}),
    ],
)

row = html.Div(
    [
    dbc.Row([
        dbc.Col(children=[title,web_des]),
        dbc.Col(
            children=[apple_plot], width=7, style={'marginLeft': 5})
        ]),
    ],
    style={'marginBottom': 50, 'marginTop': 25}
)

app.layout = html.Div(
    [
    dbc.Container(
        children=[
            row],
        fluid = True,
    ),
    dcc.Interval(
            id='interval-component',
            interval=1*1000, # in milliseconds
            n_intervals=0
    )],
    style={'marginBottom': 50, 'marginTop': 25, 'marginLeft': 10, 'marginRight': 10}
)


@app.callback(
    Output('apple_plot', 'figure'),
    [Input('gene_list', 'value')])
def update_figure(selected_gene):
    prediction_db = sqlite3.connect("obj/prediction_20220130_3d_51markers.db")
    db_list = list(prediction_db.execute("SELECT GENENAME,PREDIC from PREDDICTION WHERE ID=?", (selected_gene,)))
    fig = get_apple_plot(db_list, selected_gene, defined_colorscale)
    return fig


@app.callback(
    Output('insitu', 'src'),
    [Input('gene_list', 'value')])
def update_insitu(selected_gene):
    if( selected_gene not in insitu_gene):
        insitu_url = 'assets/insitu_NA.png'
    else:
        insitu_url = 'assets/' + selected_gene + ".png"
    return insitu_url

@app.callback(
    Output('gene_name', 'children'),
    [Input('gene_list', 'value')])
def update_insitu(selected_gene):
    try:
        gene_name = id_name_dic[selected_gene]
    except:
        gene_name = "Invalid gene name"
    return gene_name


@app.callback(
    Output('tbl', 'data'),
    [Input('gene_list', 'value')])
def update_table(selected_gene):
    try:
        data = get_similar_genes(selected_gene, cor_df, id_name_dic)
    except:
        data = None
    return data



server = app.server

if __name__ == "__main__":
    application.run(debug=True, port=8080)
