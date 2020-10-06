#!/usr/bin/env python3

import os
import logging
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

auth_url = config.conf.get("HTTP", "auth_url")
auth_port = config.conf.get("HTTP", "auth_port")
if not auth_url.startswith("http"):
    auth_url = f"http://{auth_url}"
os.environ["TOKENINFO_URL"] = f"{auth_url}:{auth_port}/tokeninfo"
basedir = os.path.abspath(os.path.dirname(__file__))
application = connexion.FlaskApp("JAWS_REST", specification_dir=basedir)
application.add_api("swagger.rest.yml")

application.app.config["SQLALCHEMY_DATABASE_URI"] = "%s://%s:%s@%s:%s/%s" % (
    config.conf.get("DB", "dialect"),
    config.conf.get("DB", "user"),
    quote_plus(config.conf.get("DB", "password")),
    config.conf.get("DB", "host"),
    config.conf.get("DB", "port"),
    config.conf.get("DB", "db"),
)
application.app.config["SQLALCHEMY_ECHO"] = False
application.app.config["SQLALCHEMY_TRACK_MODIFICATIONS"] = False
application.app.config["SQLALCHEMY_ENGINE_OPTIONS"] = {
    "pool_pre_ping": True,
    "pool_recycle": 3600,
    "pool_size": 5,
    "max_overflow": 10,
    "pool_timeout": 30
}
db.init_app(application.app)

# create tables if not exists
with application.app.app_context():
    try:
        db.create_all()
        db.session.commit()
    except Exception as error:
        db.session.rollback()
        logger.exception(f"Failed to create tables: {error}")
        raise

# init RPC clients
site_rpc_params = config.conf.get_all_sites_rpc_params()
rpc_index.rpc_index = rpc_index.RPC_Index(site_rpc_params)


if __name__ == "__main__":
    # start REST server on default port (using sock file for gunicorn)
    application.run(host="0.0.0.0", debug=False)
