import pickle, random, time, warnings
from multiprocessing import cpu_count
from typing import Callable, Iterable, Literal

import requests, tqdm
from multiprocess import pool as mpp


class TqdmParallel:
    """Class for running a function in parallel with tqdm"""

    @classmethod
    def tqdm_starmap(
        cls,
        worker_fn: Callable,
        worker_args: list[tuple],
        num_workers=cpu_count() - 1,
        cleanup_fn=lambda: None,
        processes_or_threads: Literal["processes", "threads"] = "processes",
        tqdm_class=tqdm.tqdm,
        **tqdm_kwargs,
    ):
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
            partial_results: list,
            num_expected_results: int,
        ):
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
            tqdm_kwargs["bar_format"] = "{desc}{percentage:3.0f}%|{bar:25}{r_bar}"
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
                    total=len(worker_args),
                ):
                    results.append(res)
        except IndexError as e:
            if str(e) != "pop from an empty deque":
                print(e)
                return end_actions(pool, results, len(worker_args))
            else:
                return end_actions(pool, results, len(worker_args))
        except KeyboardInterrupt:
            return end_actions(pool, results, len(worker_args))
        except Exception as e:
            print(e)
            return end_actions(pool, results, len(worker_args))

        cleanup_fn()
        print("Parallel tasks successfully finished!")
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


def request_worker(
    url: str, method: Literal["GET", "POST"] = "GET", request_kwargs: dict | None = None
) -> str | requests.HTTPError:
    if request_kwargs is None:
        request_kwargs = {}
    if "timeout" not in request_kwargs:
        request_kwargs["timeout"] = 5

    def inner_main():
        with requests.request(method, url, **request_kwargs) as r:
            if r.status_code == 200:
                return r.text
            else:
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
    num_concurrent: int | Literal["num_urls", "num_cores"] = "num_cores",
    request_kwargs: dict | None = None,
    prog_bar_kwargs: dict | None = None,
    method: Literal["GET", "POST"] = "GET",
):
    if prog_bar_kwargs is None:
        prog_bar_kwargs = {}

    if request_kwargs is None:
        request_kwargs = {}
    match num_concurrent:
        case "num_urls":
            num_concurrent = len(urls)
        case "num_cores":
            num_concurrent = cpu_count()
        case _:
            assert isinstance(
                num_concurrent, int
            ), "num_concurrent must be an int if not 'num_urls' or 'num_cores'"

    return TqdmParallel.tqdm_starmap(
        request_worker,
        [(url, method, request_kwargs) for url in urls],
        num_workers=num_concurrent,
        **prog_bar_kwargs,
    )


if __name__ == "__main__":
    # urls = [f"https://camp.bicnirrh.res.in/seqDisp.php?id=CAMPSQ{i}" for i in range(1, 24815)]
    # results = parallel_requests(urls, num_concurrent=25)
    # with open("results.pkl", "wb") as f:
    #     pickle.dump(results, f)
    results = pickle.load(open("results.pkl", "rb"))
