import jaws_central.jaws_constants


def test_run_status_msg():

    valid_run_states = [
        "created",
        "upload queued",
        "uploading",
        "upload failed",
        "upload inactive",
        "upload complete",
        "ready",
        "missing input",
        "submitted",
        "submission failed",
        "queued",
        "running",
        "succeeded",
        "failed",
        "finished",
        "aborting",
        "aborted",
        "download queued",
        "downloading",
        "download failed",
        "download inactive",
        "download complete",
    ]

    for status in valid_run_states:
        assert status in jaws_central.jaws_constants.run_status_msg
        description = jaws_central.jaws_constants.run_status_msg[status]
        assert isinstance(description, str)
        assert len(description)

    for status in jaws_central.jaws_constants.run_status_msg:
        assert status in valid_run_states
