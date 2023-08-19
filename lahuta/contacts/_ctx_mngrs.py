import contextlib
from typing import Any, Generator

import joblib


@contextlib.contextmanager
def tqdm_joblib(tqdm_object: Any) -> Generator[Any, None, None]:
    """Context manager to patch joblib to report into tqdm progress bar given as argument."""

    class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):  # type: ignore
        """Class to patch joblib to report into tqdm progress bar given as argument."""

        def __call__(self, *args: Any, **kwargs: Any) -> Any:
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)  # type: ignore

    old_batch_callback = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old_batch_callback
        tqdm_object.close()
