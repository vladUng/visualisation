
import dash
from plotlyflask import init_app
import dash_core_components as dcc
import dash_html_components as html

import webbrowser
from threading import Timer


from plotlyflask.plotlydash.main import app

from plotlyflask.plotlydash.features import gene_diff as gf
from plotlyflask.plotlydash.features import gene_viz as gv
from plotlyflask.plotlydash.features import manyfold as mf


index_page = html.Div([
    html.H1("Welcome to JBU's Data visualisation!"),
    html.H2("Below are the options available."),
    dcc.Link('Gene Viz', href='/gene-vis'),
    html.Br(),
    dcc.Link('Man ifold', href='/gene-diff'),
    html.Br(),
    dcc.Link('Manifold', href='/manyfold'),
    dcc.Loading(
        id="loading-1",
        type="default",
        children=html.Div(id="loading-output-1")
    ),
])

app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    html.Div(id='page-content'),
])


@app.callback(dash.dependencies.Output('page-content', 'children'),
              [dash.dependencies.Input('url', 'pathname')])
def display_page(pathname):
    if pathname == '/gene-vis':
        return gv.layout
    elif pathname == '/gene-diff':
        return gf.layout
    elif pathname == '/manyfold':
        return mf.layout
    else:
        return index_page
    # You could also return a 404 "URL not found" page here


def open_browser():
    # webbrowser.open_new('http://127.0.0.1:8080/dashapp')
    webbrowser.open_new('http://127.0.0.1:8080/gene-diff')

if __name__ == "__main__":
    #     Timer(1, open_browser).start()
    app.server.run(host='0.0.0.0', port=8080, debug=False)
