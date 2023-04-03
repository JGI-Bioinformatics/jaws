from jaws_site import queue_wait as slurm_queue_wait
from .conftest import MockLogger
from datetime import datetime, timedelta


def test_check_queue_wait(monkeypatch):
    def mock_sbatch_test_only():
        time1 = datetime.now() + timedelta(minutes=11)
        time1s = time1.strftime("%Y-%m-%dT%H:%M:%S")
        time2 = datetime.now() + timedelta(minutes=21)
        time2s = time2.strftime("%Y-%m-%dT%H:%M:%S")
        result = {
            "small": "",
            "medium": f"sbatch: Job 1239738 to start at {time1s} (...) nid00311 in partition genepool_shared",
            "large": "",
            "xlarge": f"sbatch: Job 3909868 to start at {time2s} (...) exvivo014 in partition exvivo",
        }
        return result

    monkeypatch.setattr(slurm_queue_wait, "_sbatch_test_only", mock_sbatch_test_only)

    logger = MockLogger()
    actual = slurm_queue_wait.check_queue_wait(logger)
    assert "small" not in actual
    assert "medium" in actual
    assert int(actual["medium"]) >= 600
    assert "large" not in actual
    assert "xlarge" in actual
    assert int(actual["xlarge"]) >= 1200
