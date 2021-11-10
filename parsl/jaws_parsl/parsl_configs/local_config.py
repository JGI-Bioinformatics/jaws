from parsl.config import Config
from parsl.providers import LocalProvider
from parsl.executors import HighThroughputExecutor
from parsl.addresses import address_by_hostname

CONFIG_LOCAL = Config(
    executors=[
        HighThroughputExecutor(
            label='local',
            address=address_by_hostname(),
            provider=LocalProvider(),
        ),
    ],
)
