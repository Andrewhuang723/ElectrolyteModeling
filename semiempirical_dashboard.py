import dash
from dash import Dash, html, dcc, Input, Output, callback, State
import pandas as pd
import plotly.express as px
import plotly.figure_factory as ff
import plotly.graph_objects as go
import configparser
import psycopg2
from sqlalchemy import create_engine
from Semiemperical import *
from query import system_chosen


app = Dash(__name__)

app.layout = html.Div(className='container',
    children=[html.Div([
        html.H1("Foxconn Battery Cell Analysis Lab"),
        html.H2("Semi-empirical Model For Gel Polymer Electrolyte"),
    ]),
    html.Div([
        html.Div([
            html.H3("Functions"),
            dcc.Dropdown(
                options=["Weito2020", "Kim2011", "Landesfeind2019"],
                value="Weito2020",
                id="function-dropdown",
                style={'width': '75%', 'color': 'black'} 
            ),
            
            html.H3("Solvent 1"),
            dcc.Dropdown(
                options=["EC", "PC", "EMC", "DMC", "PEGDM-400", "PEGDM-1000", "PEGDM-2500", "PEGDME"],
                value="EC",
                id="solvent1-dropdown",
                style={'width': '75%', 'color': 'black'} 
            ),

            html.H3("Solvent 2"),
            dcc.Dropdown(
                options=["EC", "PC", "EMC", "DMC", "PEGDM-400", "PEGDM-1000", "PEGDM-2500", "PEGDME"],
                value="PC",
                id="solvent2-dropdown",
                style={'width': '75%', 'color': 'black'} 
            ),

            html.H3("Salt"),
            dcc.Dropdown(
                options=["LiTDI", "LiAsF6", "LiTFSI", "LiBF4", "LiPF6", "TPSB-"],
                value="LiPF6",
                id="salt-dropdown",
                style={'width': '75%', 'color': 'black'} 
            ),
            ### Model run button
            html.Button(
                children='RUN', 
                id='run-model-button', 
                n_clicks=0,
                style={'background-color': 'blue', 'color': 'white'}),
            html.Div(id='output-model-file'),
        ], style={'flex': 1}
        ),

        html.Div([
            html.Div([
            html.Label("Online DB Testing results"),
            html.Table([
                html.Tr([html.Td(['Count']), html.Td(id='Testing-Count')]),
                html.Tr([html.Td(['R2']), html.Td(id='Testing-R2')]),
                html.Tr([html.Td(['MSE']), html.Td(id='Testing-MSE')]),
                html.Tr([html.Td(['Accuracy']), html.Td(id='Testing-Accuracy')]),
                html.Tr([html.Td(['Parameter']), html.Td(id='Parameter')]),
            ])],
            style={'margin-top': '10px'}),
        
            html.Div([
            html.Label("Online Testing"),
            dcc.Graph(
                id="Testing-plot",
                style={'width': '30vh', 'height': '30vh'}
            )],
        style={'flex': 1, 'text-align': 'center', 'margin-top': '10px'}
        )],
        style={'width': '50%', 'flex': 1, 'text-align': 'center'}
        )
        ,
        html.Div([
            html.H3("Concentration Profile"),
            html.H4("Temperature"),
            dcc.Input(
                id="temperature-input",
                type="number",
                value=273.15, min=250, max=400   
            ),
            dcc.Graph(
                id="concenctration-profile",
                style={'width': '80', 'height': '50vh'}
            )],
        style={'flex': 1, 'text-align': 'center', 'margin-top': '-50px'}
        )
    ],
    style={'display': 'flex'}
    ),
])


def calculate_metrics(ypred, ytest) -> dict:
    R2 = r2_score(y_pred=ypred, y_true=ytest)
    MSE = mean_squared_error(y_pred=ypred, y_true=ytest)
    error = np.abs(ypred - ytest) / ytest
    Accuracy = len(np.where(error < 0.1)[0]) / len(ytest)
    return {
        "R2": R2,
        "MSE": MSE,
        "Accuracy": Accuracy
    }


def R2_plot_plotly(ypred, ytest):
    fig = px.scatter(x=ypred, y=ytest)
    lmin = min(np.min(ypred), np.min(ytest)) - 0.01
    lmax = max(np.max(ypred), np.max(ytest)) + 0.01
    lim = [lmin, lmax]
    fig.add_trace(px.line(x=lim, y=lim).data[0])
    fig.update_xaxes(title="Conductivity test (mS/cm)", range=[lmin, lmax])
    fig.update_yaxes(title="Conductivity_pred (mS/cm)", range=[lmin, lmax])
    fig.update_layout(
        margin={'l': 40, 'b': 40, 't': 0, 'r': 0}, 
        hovermode='closest',
        width=300,
        height=300)
    return fig   


def concentration_profile_plot(salt_conc, temp, param, conductivity_pred):
    # Create line plot
    x_steps = np.linspace(0.1, 3, 50)
    fig = go.Figure()

    # Add line plot
    fig.add_trace(go.Scatter(x=x_steps, y=temp, mode='lines', name=f'T={temp}'))

    # Add scatter plot on top of the line plot
    fig.add_trace(go.Scatter(x=x_scatter, y=conductivity_pred, mode='markers', name='Scatter Plot'))


@callback(
    Output("function-dropdown", "value"),
    Input("function-dropdown", "options"),
)
def set_function(function_name):
    return function_name


@callback(
    Output("solvent1-dropdown", "value"),
    Output("solvent2-dropdown", "value"),
    Output("salt-dropdown", "value"),
    Input("solvent1-dropdown", "options"),
    Input("solvent2-dropdown", "options"),
    Input("salt-dropdown", "options"),
)
def set_system(solvent1, solvent2, salt):
    return solvent1, solvent2, salt


@callback(
    Output("Testing-plot", "figure"),
    Output("Testing-Count", "children"),
    Output("Testing-R2", "children"),
    Output("Testing-MSE", "children"),
    Output("Testing-Accuracy", "children"),
    Output("Parameter", "children"),
    Input("run-model-button", "n_clicks"),
    State("function-dropdown", "value"),
    State("solvent1-dropdown", "value"),
    State("solvent2-dropdown", "value"),
    State("salt-dropdown", "value"),
    State("output-model-file", "children"),
)
def process_function(n_clicks, function_name, solvent1, solvent2, salt, model_null):

    ## system
    df = system_chosen(s1=solvent1, s2=solvent2, salt=salt)
    
    seed = 100
    function_info = conductivity_function_chosen(function_name, seed=seed)
    curve_fit_function = function_info["curve_fit_function"]
    initial_guess = function_info["initial_guess"]

    X, Y = parse_dataset(df)
    
    try:
        xtrain, xtest, ytrain, ytest = train_test_split(X, Y, test_size=0.2, random_state=10)
    except ValueError:
        print("Not enough samples for train test split")
        
    # xtrain, xtest are tuple in the form of (solvent1, solvent2, salt, temp)

    
    parameter = run_fitting(xtrain, ytrain, 
                            curve_fit_function=curve_fit_function, 
                            initial_guess=initial_guess)
    ypred = evaluate(xtest, ytest, curve_fit_function=curve_fit_function, parameter=parameter)
    r2_fig = R2_plot_plotly(ypred=ypred, ytest=ytest)
    results = calculate_metrics(ypred=ypred, ytest=ytest)
    
    return r2_fig, len(ytest), results["R2"], results["MSE"], results["Accuracy"], parameter


@callback(
    Output("temperature-input", "value"),
    Input("temperature-input", "value")
)
def set_temperature(temperature):
    return temperature

@callback(
    Output("concenctration-profile", "figure"),
    Input("function-dropdown", "value"),
    Input("temperature-input", "value"),
    Input("Parameter", "children"),
    State("solvent1-dropdown", "value"),
    State("solvent2-dropdown", "value"),
    State("salt-dropdown", "value"),
)
def concentration_profile(function_name, temperature, parameter, solvent1, solvent2, salt):
    function_info = conductivity_function_chosen(function_name)
    function = function_info["function"]

    # Create line plot
    fig = go.Figure()

    xs = np.linspace(0, 3, 50)
    ys = [function(0, 0, x, temperature, *parameter) for x in xs]
    # Add line plot
    fig.add_trace(go.Scatter(x=xs, y=ys, mode='lines', name='Concentration Profile'))

    # Add scatter plot on top of the line plot
    # fig.add_trace(go.Scatter(x=x_scatter, y=y_scatter, mode='markers', name='Scatter Plot'))

    # Customize the layout
    fig.update_layout(title=f'Concentration profile of {solvent1}/{solvent2}/{salt} @ Temp = {temperature}',
                    xaxis_title=f'{salt} (Mol / Kg Polymer)',
                    yaxis_title='Conductivtiy (mS/cm)',
                    width=500, height=300
                    )
    return fig


if __name__ == "__main__":
    app.run_server(debug=True, port=8050)