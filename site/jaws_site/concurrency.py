import functools
from collections.abc import AsyncIterator, Callable, Iterator
from typing import ParamSpec, TypeVar

import anyio

T = TypeVar("T")
P = ParamSpec("P")


async def run_in_threadpool(
    func: Callable[P, T], *args: P.args, **kwargs: P.kwargs
) -> T:
    """Run a function in a separate thread.

    Many of our processes are still synchronous so to stop them from
    interfering with the main event loop we will run them in
    background threads.
    """
    if kwargs:
        func = functools.partial(func, **kwargs)
    return await anyio.to_thread.run_sync(func, *args)


class _StopIteration(Exception):
    pass


def _next(iterator: Iterator[T]) -> T:  # type: ignore
    """Coerce threadpool iterator into different exception type.

    We cannot raise `StopIteration` from within the threadpool
    iterator and catch it outside of that context, so we coerce a
    different exception.
    """
    try:
        return next(iterator)
    except StopIteration:
        pass


async def iterate_in_threadpool(iterator: Iterator[T]) -> AsyncIterator[T]:
    while True:
        try:
            yield await anyio.to_thread.run_sync(_next, iterator)
        except _StopIteration:
            break
