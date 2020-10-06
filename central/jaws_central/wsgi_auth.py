#!/usr/bin/env python3

"""JAWS Central CLI"""

import os
import logging
import connexion
from urllib.parse import quote_plus
from jaws_central import config, log
from jaws_central.models_fsa import db

JAWS_LOG_ENV = "JAWS_CENTRAL_LOG"
JAWS_CONFIG_ENV = "JAWS_CENTRAL_CONFIG"


if __name__ == "__main__":
    if JAWS_LOG_ENV not in os.environ:
        raise SystemExit(f"Required env undefined: {JAWS_LOG_ENV}")
    if JAWS_CONFIG_ENV not in os.environ:
        raise SystemExit(f"Required env undefined: {JAWS_CONFIG_ENV}")

    log_file = os.environ[JAWS_LOG_ENV]
    logger = log.setup_logger(__package__, log_file, log_level)
    config_file = os.environ[JAWS_CONFIG_ENV]
    conf = config.Configuration(config_file)
    if conf:
        logger.debug(f"Config using {config_file}")

    """Start JAWS OAuth server"""
    logger = logging.getLogger(__package__)
    logger.debug("Initializing OAuth server")
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

    # define port
    port = int(config.conf.get("HTTP", "auth_port"))  # defaults to 3000

    # start OAuth server
    logger.debug(f"Starting OAuth server on port {port}")
    #connex.run(host="0.0.0.0", port=port, debug=False)
    connex.run()
