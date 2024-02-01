import logging

from jaws_common import messages

from jaws_site.task_logger import TaskLogger


def save_task_log(params, **kwargs):
    config = kwargs.get("config")
    session = kwargs.get("session")
    logger = kwargs.get("logger")
    task_logger = TaskLogger(config, session, logger=logger)
    task_logger.save(params)


# dispatch table
operations = {
    "task_logger": {
        "function": save_task_log,
        "required_params": ["status", "cromwell_run_id"],
    },
}


class Consumer:
    def __init__(self, config, session, **kwargs):
        self.config = config
        self.session = session
        self.logger = kwargs.get("logger", None)
        if self.logger is None:
            self.logger = logging.getLogger()
        site_id = self.config.get("SITE", "id")
        deployment = self.config.get("SITE", "deployment")
        self.queue = f"jaws_{site_id}_{deployment}"
        rmq_config = self.config.get_section("RMQ")
        self.consumer = messages.Consumer(
            config=rmq_config,
            session=self.session,
            logger=self.logger,
            operations=operations,
            queue=self.queue,
        )

    def consume(self):
        self.consumer.consume()
