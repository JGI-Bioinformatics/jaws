from collections.abc import Callable, Sequence
from typing import Any, ParamSpec

from .concurrency import run_in_threadpool

P = ParamSpec("P")


class BackgroundTask:
    """Run function in a background thread.

    A background task is a function that takes another function and
    its inputs and runs them in a background thread in the running
    threadpool. Python is not an asynchronous language by default so a
    long running function (one that does IO synchronously or lots of
    computation) has the potential to block the rest of the program,
    making it counterproductive in an event loop.
    """
    def __init__(
            self, func: Callable[P, Any], *args: P.args, **kwargs: P.kwargs
    ) -> None:
        self.func = func
        self.args = args
        self.kwargs = kwargs

    async def __call__(self) -> None:
        await run_in_threadpool(self.func, *self.args, **self.kwargs)


class BackgroundTasks(BackgroundTask):
    def __init__(self, tasks: Sequence[BackgroundTask] | None) -> None:
        self.tasks = list(tasks) if tasks else []

    def add_task(
        self, func: Callable[P, Any], *args: P.args, **kwargs: P.kwargs
    ) -> None:
        """Add function to be executed in threadpool.

        A function passed to this method is offloaded to a background
        thread and executed independently of the main program's event
        loop. This prevents potentially blocking code from holding up
        asynchronous tasks.
        """
        task = BackgroundTask(func, *args, **kwargs)
        self.tasks.append(task)

    async def __call__(self) -> None:
        for task in self.tasks:
            await task()
