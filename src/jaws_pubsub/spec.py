from abc import ABC, abstractmethod
from collections.abc import Callable
from dataclasses import dataclass
from enum import Enum
from logging import Logger
from typing import Any, TypeVar

DecoratedCallable = TypeVar("DecoratedCallable", bound=Callable[..., Any])

class Priority(Enum):
    LOW = 1
    NORMAL = 5
    HIGH = 10


@dataclass
class Headers:
    priority: Priority
    task_type: str


class MessageClient(ABC):
    """Instantiates and object for message passing with external broker.

    :param configuration: Object containing auth parameters for broker
    :type configuration: dict[str,str]
    :param logger: Logging Transport
    :type logger: class:`logging.Logger`"""

    def __init__(
        self,
        configuration: dict[str,str],
        logger: Logger,
    ) -> None:
        """Create a message client and connect to a broker.

        :param configuration: Object containing parameters for broker
        :type configuration: dict[str,str]

        :param logger: Logging Transport
        :type logger: class:`logging.Logger`
        """
        raise NotImplementedError

    @abstractmethod
    def _connect(self) -> None:
        raise NotImplementedError

    @abstractmethod
    def _init_connection_parameters(self) -> None:
        raise NotImplementedError

    @abstractmethod
    def check_connection(self) -> None:
        raise NotImplementedError

    @abstractmethod
    def close(self) -> None:
        raise NotImplementedError

    @abstractmethod
    def declare_queue(
        self,
        queue_name: str,
        exclusive: bool,
        max_priority: int,
    ) -> None:
        """Declare queue, set its exclusivity, and give it priority.

        By default the `queue_name` is '', which means naming of the
        queue is left up to this method as opposed to the
        caller. `max_priority` has three levels: LOW(1), NORMAL(5),
        and HIGH(10). By default `max_priority` is set to NORMAL(5),
        meaning messages in this queue cannot have higher priority
        than normal.

        :param queue_name: The name of the queue 
        :type queue_name: str
        :param exclusive: Whether the queue will be exclusive to a receiver or not
        :type exclusive: bool
        :param max_priority: Priority given to messages in the queue 
        :type max_priority: int

        :raises NotImplementedError: if this method does not handle
            `queue_name` being an empty string.
        """
        raise NotImplementedError

    @abstractmethod
    def declare_exchange(self, exchange_name: str, exchange_type: str) -> None:
        """Declare the exchange and its type.

        :param exchange_name: The name of the exchange
        :type exchange_name: str
        :param exchange_type: Type of exchange
        :type exchange_type: str

        :raises NotImplementedError: if exchange name or type is not defined"""
        raise NotImplementedError

    @abstractmethod
    def bind_queue(self, exchange_name: str, queue_name: str, routing_key: str) -> None:
        """Bind queue to an exchange."""
        raise NotImplementedError

    @abstractmethod
    def unbind_queue(
        self, exchange_name: str, queue_name: str, routing_key: str
    ) -> None:
        """Unbind a queue and exchange."""
        raise NotImplementedError


class MessageSender(MessageClient):
    @abstractmethod
    def encode_message(self, body: object, encoding_type: str) -> str:
        """Encodes message into format recognized by the exchange.
        :param body: object to be encoded for the exchange 
        :type body: object
        :param encoding_type: the type of the encoding
        :type encoding_type: str

        :raises NotImplementedError: if declared encoding_type is not implemented 
        """
        raise NotImplementedError

    @abstractmethod
    def send_message(
        self,
        exchange_name: str,
        routing_key: str,
        body: object,
        headers: Headers | None,
    ) -> None:
        raise NotImplementedError


class MessageReceiver(MessageClient):

    @abstractmethod
    def decode_message(self, body: bytes | str) -> dict[str, str]:
        """Decodes message from exchange into a python dictionary.

        :param body: The body of the message received from the broker
        :type body: bytes, str

        :raises NotImplementedError: If the body is in an unrecognized format"""
        raise NotImplementedError

    @abstractmethod
    def get_message(self, queue_name: str, auto_ack: bool) -> dict[str, str] | None:
        """Returns message from queue as a python object."""
        raise NotImplementedError

    @abstractmethod
    def consume_messages(self, queue: str, callback: Callable[..., None]) -> None:
        """Gathers messages from the queue and applies a callback function to them."""
        raise NotImplementedError

    @abstractmethod
    def cancel_consumer(self) -> None:
        """Disconnects consumer from the message queue."""
        raise NotImplementedError
