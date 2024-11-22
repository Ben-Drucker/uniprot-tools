"""Tools to do things in parallel with progress bars"""

import gzip, io, os, pickle, random, requests, time, tqdm, traceback, warnings
from multiprocess import pool as mpp
from multiprocessing import cpu_count
from termcolor import colored
from typing import Callable, Generator, Iterable, Literal, TypeVar

T = TypeVar("T")


class TqdmParallel:
    """Class for running a function in parallel with tqdm"""

    @classmethod
    def tqdm_starmap(
        cls,
        worker_fn: Callable[..., T],
        worker_args: list[tuple] | Generator,
        num_workers=cpu_count() - 1,
        cleanup_fn=lambda: None,
        processes_or_threads: Literal["processes", "threads"] = "processes",
        tqdm_class=tqdm.tqdm,
        n_args: int | None = None,
        **tqdm_kwargs,
    ) -> list[T] | None:
        """Run a function in parallel with tqdm

        Parameters
        ----------
        ``worker_fn`` :
            The function to run in parallel
        ``worker_args`` :
            List of tuples of arguments to pass to the function. Each tuple is passed to the \
                ``worker_fn``'s ``*args``.
        ``num_workers`` :
            The number of workers to use. Defaults to the number of cores - 1.
        ``cleanup_fn`` :
            The function to call when the workers are finished. Defaults to a no-op function.
        ``processes_or_threads`` :
            Whether to use multiprocessing or threading. Defaults to multiprocessing.
        ``tqdm_class`` :
            The tqdm class to use. Defaults to ``tqdm.tqdm``, but could also be, for example, \
                ``tqdm.tqdm_notebook``.
        ``tqdm_kwargs`` :
            Additional keyword arguments to pass to the tqdm class. Some possible options include::

                desc=
                leave=
                ncols=
                ascii=
                disable=
                unit=
                unit_scale=
                dynamic_ncols=
                smoothing=
                position=

        """

        def end_actions(
            pool: mpp.Pool | mpp.ThreadPool | None,
            partial_results: list[T],
            num_expected_results: int,
        ) -> list[T] | None:
            tqdm.tqdm()
            if pool is not None:
                pool.terminate()
                pool.close()
                cleanup_fn()
                print(
                    f"Returning partial results ([{len(partial_results)}/{num_expected_results}])"
                )
                return partial_results  # , worker_args[len(partial_results):]

        results = []
        pool = None
        if "smoothing" not in tqdm_kwargs:
            tqdm_kwargs["smoothing"] = 0.05
        if "leave" not in tqdm_kwargs:
            tqdm_kwargs["leave"] = True
        if "desc" not in tqdm_kwargs:
            tqdm_kwargs["desc"] = "Progress: "
        if "bar_format" not in tqdm_kwargs:
            tqdm_kwargs["bar_format"] = "{desc} {percentage:3.0f}% {bar:25}"

        if isinstance(worker_args, list):
            n_args = len(worker_args)
        else:
            if isinstance(worker_args, Generator):
                if n_args is None:
                    raise ValueError(
                        "If `worker_args` is a Generator, the number of args (`n_args`) must be"
                        " provided to the current function, `tqdm_starmap`."
                    )

        try:
            match processes_or_threads:
                case "processes":
                    cls_ = mpp.Pool
                case "threads":
                    cls_ = mpp.ThreadPool
            with cls_(num_workers) as pool:
                for res in tqdm_class(
                    cls._istarmap(pool, worker_fn, worker_args),
                    **tqdm_kwargs,
                    total=n_args,
                ):
                    results.append(res)
        except IndexError as e:
            if str(e) != "pop from an empty deque":
                print(e)
                res = end_actions(pool, results, n_args)
            else:
                return end_actions(pool, results, n_args)
        except KeyboardInterrupt:
            return end_actions(pool, results, n_args)
        except Exception as e:
            print(colored("\n\t".join(traceback.format_exc().split("\n")), color="red"))
            return end_actions(pool, results, n_args)

        cleanup_fn()
        # print("Parallel tasks successfully finished!")
        return results

    @classmethod
    def _istarmap(
        cls, p: mpp.Pool | mpp.ThreadPool, func: Callable, iterable: Iterable, chunksize=1
    ):
        """
        Notes
        -----
        Copied from from https://stackoverflow.com/a/57364423"""
        p._check_running()  # type: ignore
        if chunksize < 1:
            raise ValueError("Chunksize must be 1+, not {0:n}".format(chunksize))

        task_batches = mpp.Pool._get_tasks(func, iterable, chunksize)  # type: ignore
        result = mpp.IMapIterator(p)
        p._taskqueue.put(  # type: ignore
            (
                p._guarded_task_generation(result._job, mpp.starmapstar, task_batches),  # type: ignore
                result._set_length,  # type: ignore
            )
        )
        return (item for chunk in result for item in chunk)


def _request_worker(
    url: str,
    method: Literal["GET", "POST"] = "GET",
    post_data: dict | None = None,
    stream_and_dump_path: str | None = None,
    request_kwargs: dict | None = None,
) -> str | requests.HTTPError:
    if request_kwargs is None:
        request_kwargs = {}
    if "timeout" not in request_kwargs:
        request_kwargs["timeout"] = 5

    assert method in ["GET", "POST"], f"method must be GET or POST. It was {method}."

    def inner_main():
        if method == "POST":
            if "data" in request_kwargs and post_data is not None:
                raise ValueError("Cannot use both `post_data` and `request_kwargs['data']`")

            request_kwargs["data"] = post_data
        with requests.request(method, url, **request_kwargs) as r:
            try:
                r.raise_for_status()
                if stream_and_dump_path is None:
                    # return either text or bytes
                    try:
                        with gzip.GzipFile(fileobj=io.BytesIO(r.content)) as f:
                            decoded = f.read().decode("utf-8")
                            return decoded
                    except gzip.BadGzipFile as e:
                        if "not a gzipped file" not in str(e).lower():
                            raise e from None
                        return r.text
                else:
                    with open(stream_and_dump_path, "wb") as f:
                        for chunk in r.iter_content(chunk_size=2**10):
                            f.write(chunk)
                    return stream_and_dump_path
            except requests.exceptions.HTTPError:
                msg = f"Got status code {r.status_code} from {url}. Response: {r.text}"
                warnings.warn(msg)
                return requests.HTTPError(msg)

    while True:
        try:
            return inner_main()
        except requests.exceptions.ReadTimeout:
            time.sleep(10 * random.random())
            continue
        except requests.exceptions.ConnectionError:
            time.sleep(10 * random.random())
            continue


def parallel_requests(
    urls: list[str],
    post_datas: list[dict | None] | None = None,
    num_concurrent: int | Literal["num_urls", "num_cores"] = "num_cores",
    request_kwargs: dict | None = None,
    prog_bar_kwargs: dict | None = None,
    method: Literal["GET", "POST"] = "GET",
    stream_and_dump_dir: str | None = None,
) -> list[str | requests.HTTPError] | None:
    """Do HTTP(S) requests in parallel for GET or POST

    Parameters
    ----------
    ``urls`` :
        List of URLs to request
    ``post_datas`` :
        If ``method`` is POST, list of dictionaries of POST data
    ``num_concurrent`` :
        The number of concurrent requests. If 'num_urls', the number of requests is equal to \
            the number of URLs. If 'num_cores', the number of requests is equal to the number of cores on the machine.
    ``request_kwargs`` :
        Keyword arguments to pass to requests.request
    ``prog_bar_kwargs`` :
        Keyword arguments to pass to tqdm.tqdm (the progress bar)
    ``method`` :
        Either 'GET' or 'POST'
    ``stream_and_dump_dir`` :
        If not None, stream the response and save to this directory.

    Returns
    -------
        If `None`, then there was an error. Otherwise, a list of strings containing response text \
            or requests.HTTPError objects. If ``stream_and_dump_dir`` is not None, \
            this variable will simply be returned, not response text.
    """
    if prog_bar_kwargs is None:
        prog_bar_kwargs = {}

    if stream_and_dump_dir is not None:
        if not os.path.exists(stream_and_dump_dir):
            os.makedirs(stream_and_dump_dir)

    if request_kwargs is None:
        request_kwargs = {}
    request_kwargs.update({"stream": True})
    match num_concurrent:
        case "num_urls":
            num_concurrent = len(urls)
        case "num_cores":
            num_concurrent = cpu_count()
        case _:
            assert isinstance(
                num_concurrent, int
            ), "num_concurrent must be an int if not 'num_urls' or 'num_cores'"

    if stream_and_dump_dir is not None:
        paths = [
            os.path.join(stream_and_dump_dir, f"download-{i}.download") for i in range(len(urls))
        ]
    else:
        paths = [None] * len(urls)

    if post_datas is None:
        zip_iter = zip(urls, [None] * len(urls), paths)
    else:
        zip_iter = zip(urls, post_datas, paths)

    return TqdmParallel.tqdm_starmap(
        _request_worker,
        [
            (url, method, post_data, dump_path, request_kwargs)
            for url, post_data, dump_path in zip_iter
        ],
        num_workers=num_concurrent,
        processes_or_threads="threads",
        **prog_bar_kwargs,
    )
