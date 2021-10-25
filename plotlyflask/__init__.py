"""Initialize Flask app."""
from flask import Flask

# from visualisation.plotlyflask.plotlydash.dashboard import init_dashboard


def init_app():
    """Construct core Flask application."""
    app = Flask(__name__, instance_relative_config=True)
    app.config.from_object('config.Config')

    with app.app_context():
        # Import parts of our core Flask app
        from . import routes

        # Import Dash application
        # from .plotlydash.dashboard import init_dashboard
        from .plotlydash.main import app

        # app = init_dashboard(app)

        return app
