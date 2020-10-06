#!/usr/bin/env python3

import os
import connexion
from urllib.parse import quote_plus
from jaws_central import config, log
from jaws_central.models_fsa import db

JAWS_LOG_ENV = "JAWS_CENTRAL_LOG"
JAWS_LOG_LEVEL_ENV = "JAWS_CENTRAL_LOG_LEVEL"
JAWS_CONFIG_ENV = "JAWS_CENTRAL_CONFIG"

if JAWS_LOG_ENV not in os.environ:
    raise SystemExit(f"Required env undefined: {JAWS_LOG_ENV}")
if JAWS_LOG_LEVEL_ENV not in os.environ:
    raise SystemExit(f"Required env undefined: {JAWS_LOG_LEVEL_ENV}")
if JAWS_CONFIG_ENV not in os.environ:
    raise SystemExit(f"Required env undefined: {JAWS_CONFIG_ENV}")

log_file = os.environ[JAWS_LOG_ENV]
log_level = os.environ[JAWS_LOG_LEVEL_ENV]
config_file = os.environ[JAWS_CONFIG_ENV]

logger = log.setup_logger(__package__, log_file, log_level)
conf = config.Configuration(config_file)
logger.debug(f"Config using {config_file}")

basedir = os.path.abspath(os.path.dirname(__file__))
connex = connexion.FlaskApp(__name__, specification_dir=basedir)
connex.add_api("swagger.auth.yml")

connex.app.config["SQLALCHEMY_DATABASE_URI"] = "%s://%s:%s@%s:%s/%s" % (
    config.conf.get("DB", "dialect"),
    config.conf.get("DB", "user"),
    quote_plus(config.conf.get("DB", "password")),
    config.conf.get("DB", "host"),
    config.conf.get("DB", "port"),
    config.conf.get("DB", "db"),
)
connex.app.config["SQLALCHEMY_ECHO"] = False
connex.app.config["SQLALCHEMY_TRACK_MODIFICATIONS"] = False
connex.app.config["SQLALCHEMY_ENGINE_OPTIONS"] = {
    "pool_pre_ping": True,
    "pool_recycle": 3600,
    "pool_size": 5,
    "max_overflow": 10,
    "pool_timeout": 30,
}
db.init_app(connex.app)

# create tables if not exists
with connex.app.app_context():
    try:
        db.create_all()
        db.session.commit()
    except Exception as error:
        db.session.rollback()
        logger.exception(f"Failed to create tables: {error}")
        raise

if __name__ == "__main__":
    # start on default port (using sock file with gunicorn)
    connex.run(host="0.0.0.0", debug=False)
