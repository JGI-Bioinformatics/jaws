from collections.abc import Callable

from jaws_pubsub.log import get_logger

from .queue import PikaMessageReceiver, PikaMessageSender
from .spec import DecoratedCallable, MessageReceiver, MessageSender

logging = get_logger("pubsub")


class JawsPubsub:
    def __init__(self, configuration: dict[str, str], logger: logging.Logger) -> None:
        self.configuration = configuration
        self.logger = logger

    def consume(
        self,
        *,
        queue_name: str = "",
        exchange_type: str = "",
        exchange_name: str = "",
        tag: str | None = None,
        exclusive: bool = False,
        max_priority: int = 5,
        consumer: MessageReceiver,
    ) -> Callable[[DecoratedCallable], DecoratedCallable]:
        """"""
        if not consumer:
            consumer = PikaMessageReceiver(tag, self.configuration, self.logger)
            consumer.declare_queue(
                queue_name=queue_name, exclusive=exclusive, max_priority=max_priority
            )
            consumer.declare_exchange(
                exchange_name=exchange_name, exchange_type=exchange_type
            )

        def decorator(
            func: DecoratedCallable, consumer: MessageReceiver = consumer
        ) -> DecoratedCallable:
            return func

        return decorator

    def publish(
        self,
        *,
        queue_name: str = "",
        exchange_type: str = "",
        exchange_name: str = "",
        producer: MessageSender | None = None,
        tag: str | None = None,
        exclusive: bool = False,
        max_priority: int = 5,
    ) -> None:
        if not producer:
            producer = PikaMessageSender(self.configuration, self.logger)
