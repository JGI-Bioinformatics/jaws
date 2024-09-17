import asyncio
import json
import logging
import ssl
import time
from collections.abc import Callable
from dataclasses import asdict
from functools import wraps
from typing import Any

import pika
from Enum import Enum
from pika.adapters.blocking_connection import BlockingChannel
from pika.exceptions import AMQPConnectionError
from pika.spec import Basic, BasicProperties

from .spec import Headers, MessageClient, MessageReceiver, MessageSender


def sync(f: Callable[..., Any]) -> Callable[..., Any]:
    @wraps(f)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        return asyncio.get_event_loop().run_until_complete(f(*args, **kwargs))

    return wrapper


class PikaMessageClient(MessageClient):
    def __init__(
        self,
        configuration: dict[str, str],
        logger: logging.Logger = logging.getLogger(__name__),
    ) -> None:
        self.logger = logger
        self.username = configuration["user"]
        self.host = configuration["host"]
        self.password = configuration["password"]
        self.port = configuration["port"]
        self.vhost = configuration["vhost"]
        self.protocol = configuration["amqp"]

        self._init_connection_parameters()
        self._connect()

    def _connect(self) -> None:
        tries = 0
        while True:
            try:
                self.connection = pika.BlockingConnection(self.parameters)
                self.channel = self.connection.channel()
                if self.connection.is_open:
                    break
            except AMQPConnectionError as e:
                time.sleep(5)
                tries += 1
                if tries == 20:
                    raise AMQPConnectionError(e)

    def _init_connection_parameters(self) -> None:
        self.credentials = pika.PlainCredentials(self.username, self.password)
        self.parameters = pika.ConnectionParameters(
            self.host,
            int(self.port),
            "/",
            self.credentials,
        )
        if self.protocol == "amqps":
            ssl_context = ssl.SSLContext(ssl.PROTOCOL_TLSv1_2)
            ssl_context.set_ciphers("ECDHE+AESGCM:!ECDSA")
            self.parameters.ssl_options = pika.SSLOptions(context=ssl_context)

    def check_connection(self) -> None:
        if not self.connection or self.connection.is_closed:
            self._connect()

    def close(self) -> None:
        self.channel.close()
        self.connection.close()

    def declare_queue(
        self, queue_name: str = "", exclusive: bool = False, max_priority: int = 5
    ) -> None:
        self.check_connection()
        self.logger.debug(f"Trying to declare queue({queue_name})...")
        self.channel.queue_declare(
            queue=queue_name,
            exclusive=exclusive,
            durable=True,
            arguments={"x-max-priority": max_priority},
        )

    def declare_exchange(self, exchange_name: str, exchange_type: Enum) -> None:
        self.check_connection()
        self.channel.exchange_declare(
            exchange=exchange_name, exchange_type=exchange_type
        )

    def bind_queue(self, exchange_name: str, queue_name: str, routing_key: str) -> None:
        self.check_connection()
        self.channel.queue_bind(
            exchange=exchange_name, queue=queue_name, routing_key=routing_key
        )

    def unbind_queue(
        self, exchange_name: str, queue_name: str, routing_key: str
    ) -> None:
        self.channel.queue_unbind(
            queue=queue_name, exchange=exchange_name, routing_key=routing_key
        )


class PikaMessageSender(PikaMessageClient, MessageSender):
    def encode_message(self, body: object, encoding_type: str = "json") -> str:
        return str(json.dumps(body))

    def send_message(
        self,
        exchange_name: str,
        routing_key: str,
        body: object,
        headers: Headers | None,
    ) -> None:
        body = self.encode_message(body=body)
        self.channel.basic_publish(
            exchange=exchange_name,
            routing_key=routing_key,
            body=body,
            properties=pika.BasicProperties(
                delivery_mode=pika.spec.PERSISTENT_DELIVERY_MODE,
                priority=headers.priority.value if headers else None,
                headers=asdict(headers) if headers else None,
            ),
        )
        self.logger.debug(
            f"Sent message. Exchange: {exchange_name}, Body: {body[:128]}"
        )


class PikaMessageReceiver(PikaMessageClient, MessageReceiver):
    def __init__(
        self, tag: str | None, configuration: dict[str, str], logger: logging.Logger
    ) -> None:
        super().__init__(configuration, logger)
        self.channel_tag = tag

    def decode_message(self, body: str | bytes) -> dict[str, str]:
        return dict(json.loads(body))

    def get_message(
        self, queue_name: str, auto_ack: bool = False
    ) -> tuple[Basic.GetOk | None, BasicProperties | None, str | None]:
        method_frame, header_frame, body = self.channel.basic_get(
            queue=queue_name, auto_ack=auto_ack
        )
        if not method_frame:
            self.logger.debug("No message returned.")
            return None, None, None
        self.logger.debug(f"{method_frame}, {header_frame}, {body}")
        return (
            (method_frame if method_frame else None),
            (header_frame if header_frame else None),
            (body if body else None),
        )

    def consume_messages(
        self,
        queue: str,
        callback: Callable[
            [BlockingChannel, Basic.Deliver, BasicProperties, bytes],
            None,
        ],
    ) -> None:
        self.check_connection()
        self.channel_tag = self.channel.basic_consume(
            queue=queue, on_message_callback=callback, auto_ack=True
        )
        self.logger.debug(" [*] Waiting for messages.")
        self.channel.start_consuming()

    def cancel_consumer(self) -> None:
        if self.channel_tag is not None:
            self.channel.basic_cancel(self.channel_tag)
            self.channel_tag = None
        else:
            self.logger.error("Do not cancel a non-existing job")


class PikaConsumer(PikaMessageReceiver):
    @sync
    async def consume(
        self,
        channel: BlockingChannel,
        method: Basic.Deliver,
        properties: BasicProperties,
        body: bytes,
    ) -> None:
        self.decode_message(body=body)
