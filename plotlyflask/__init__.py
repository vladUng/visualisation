"""Initialize Flask app."""
from flask import Flask

def init_app():
    """Construct core Flask application."""
    app = Flask(__name__, instance_relative_config=True)
    app.config.from_object('config.Config')

    with app.app_context():

        # Import Dash application
        from .plotlydash.main import app

        return app
