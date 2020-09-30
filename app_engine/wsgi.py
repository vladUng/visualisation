
from plotlyflask import init_app

import webbrowser
from threading import Timer

app = init_app()

def open_browser():
      webbrowser.open_new('http://127.0.0.1:8080/dashapp')
      
if __name__ == "__main__":
    # Timer(1, open_browser).start()
    app.run(host='0.0.0.0', port=8080, debug=True)