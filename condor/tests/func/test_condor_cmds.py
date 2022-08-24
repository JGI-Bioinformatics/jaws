from jaws_condor.htcondor_cmds import parse_condor_outputs, process_condor_q


def test_call_condor_command():

    test_data = [
        ("""38597 117760 32 834.0 37001.0 2 1 1 slot1_1@nid00333 1657822529 1657822525

""",
         [['38597', '117760', '32', '834.0', '37001.0', '2', '1', '1', 'slot1_1@nid00333', '1657822529', '1657822525']]  # noqa
         ),
         ("""

         """,  # noqa
        []  # noqa
        ),
        (
            """
            38611 2048 2 0.0 0.0 2 1 1 slot1_1@nid00333 1657839752 1657839749
            38612 2048 2 0.0 0.0 2 1 1 slot1_1@nid00362 1657839753 1657839750
            38613 2048 2 0.0 0.0 2 1 1 slot1_1@nid00414 1657839752 1657839751
            38614 2048 2 0.0 0.0 2 1 1 slot1_1@nid00301 1657839752 1657839752
            38615 2048 2 0.0 0.0 2 1 1 slot1_2@nid00333 1657839758 1657839753
            38616 2048 2 0.0 0.0 2 1 1 slot1_2@nid00362 1657839758 1657839754
            38617 2048 2 0.0 0.0 2 1 1 slot1_2@nid00414 1657839758 1657839755
            38618 2048 2 0.0 0.0 2 1 1 slot1_2@nid00301 1657839758 1657839756
            38619 2048 2 0.0 0.0 2 1 1 slot1_3@nid00333 1657839758 1657839757
            38620 2048 2 0.0 0.0 2 1 1 slot1_3@nid00362 1657839758 1657839758
            38621 2048 2 0.0 0.0 2 1 1 slot1_4@nid00333 1657839762 1657839759
            38622 2048 2 0.0 0.0 2 1 1 slot1_4@nid00362 1657839762 1657839760
            38623 2048 2 0.0 0.0 2 1 1 slot1_3@nid00414 1657839762 1657839761
            38624 2048 2 0.0 0.0 2 1 1 slot1_3@nid00301 1657839762 1657839762
            38625 2048 2 0.0 0.0 2 1 1 slot1_5@nid00333 1657839767 1657839763
            38626 2048 2 0.0 0.0 2 1 1 slot1_5@nid00362 1657839767 1657839764
            38627 2048 2 0.0 0.0 2 1 1 slot1_4@nid00414 1657839767 1657839765
            38628 2048 2 0.0 0.0 2 1 1 slot1_4@nid00301 1657839767 1657839766
            38629 2048 2 0.0 0.0 2 1 1 slot1_5@nid00301 1657839767 1657839767
            38630 2048 2 0.0 0.0 2 1 1 slot1_6@nid00333 1657839772 1657839768

            """,  # noqa
                [  # noqa
                    ['38611', '2048', '2', '0.0', '0.0', '2', '1', '1', 'slot1_1@nid00333', '1657839752', '1657839749'],
                    ['38612', '2048', '2', '0.0', '0.0', '2', '1', '1', 'slot1_1@nid00362', '1657839753', '1657839750'],
                    ['38613', '2048', '2', '0.0', '0.0', '2', '1', '1', 'slot1_1@nid00414', '1657839752', '1657839751'],
                    ['38614', '2048', '2', '0.0', '0.0', '2', '1', '1', 'slot1_1@nid00301', '1657839752', '1657839752'],
                    ['38615', '2048', '2', '0.0', '0.0', '2', '1', '1', 'slot1_2@nid00333', '1657839758', '1657839753'],
                    ['38616', '2048', '2', '0.0', '0.0', '2', '1', '1', 'slot1_2@nid00362', '1657839758', '1657839754'],
                    ['38617', '2048', '2', '0.0', '0.0', '2', '1', '1', 'slot1_2@nid00414', '1657839758', '1657839755'],
                    ['38618', '2048', '2', '0.0', '0.0', '2', '1', '1', 'slot1_2@nid00301', '1657839758', '1657839756'],
                    ['38619', '2048', '2', '0.0', '0.0', '2', '1', '1', 'slot1_3@nid00333', '1657839758', '1657839757'],
                    ['38620', '2048', '2', '0.0', '0.0', '2', '1', '1', 'slot1_3@nid00362', '1657839758', '1657839758'],
                    ['38621', '2048', '2', '0.0', '0.0', '2', '1', '1', 'slot1_4@nid00333', '1657839762', '1657839759'],
                    ['38622', '2048', '2', '0.0', '0.0', '2', '1', '1', 'slot1_4@nid00362', '1657839762', '1657839760'],
                    ['38623', '2048', '2', '0.0', '0.0', '2', '1', '1', 'slot1_3@nid00414', '1657839762', '1657839761'],
                    ['38624', '2048', '2', '0.0', '0.0', '2', '1', '1', 'slot1_3@nid00301', '1657839762', '1657839762'],
                    ['38625', '2048', '2', '0.0', '0.0', '2', '1', '1', 'slot1_5@nid00333', '1657839767', '1657839763'],
                    ['38626', '2048', '2', '0.0', '0.0', '2', '1', '1', 'slot1_5@nid00362', '1657839767', '1657839764'],
                    ['38627', '2048', '2', '0.0', '0.0', '2', '1', '1', 'slot1_4@nid00414', '1657839767', '1657839765'],
                    ['38628', '2048', '2', '0.0', '0.0', '2', '1', '1', 'slot1_4@nid00301', '1657839767', '1657839766'],
                    ['38629', '2048', '2', '0.0', '0.0', '2', '1', '1', 'slot1_5@nid00301', '1657839767', '1657839767'],
                    ['38630', '2048', '2', '0.0', '0.0', '2', '1', '1', 'slot1_6@nid00333', '1657839772', '1657839768']
                ]  # noqa
        )
        ]

    for q_output, expected in test_data:
        result = parse_condor_outputs(q_output)
        assert result == expected


def test_process_condor_q():
    columns = "ClusterId RequestMemory RequestCpus \
            CumulativeRemoteSysCpu CumulativeRemoteUserCpu \
            JobStatus NumShadowStarts JobRunCount RemoteHost JobStartDate QDate".split()

    test_data = [
        ([], [{k: 0 for k in columns}]),
        ([  # noqa
                ['38611', '2048', '2', '0.0', '0.0', '2', '1', '1', 'slot1_1@nid00333', '1657839752', '1657839749'],
                ['38612', '2048', '2', '0.0', '0.0', '2', '1', '1', 'slot1_1@nid00362', '1657839753', '1657839750'],
                ['39649', '117760', '32', '0.0', '0.0', '1', 'undefined',
                    'undefined', 'undefined', 'undefined', '1658136630'],
                ['39650', '10240', '1', '0.0', '0.0', '1', 'undefined',
                    'undefined', 'undefined', 'undefined', '1658145932'],
                ['39651', '117760', '32', '0.0', '0.0', '1', 'undefined',
                    'undefined', 'undefined', 'undefined', '1658152224'],
            ],  # noqa
            [
                {
                    'ClusterId': '38611',
                    'RequestMemory': 2.0,
                    'RequestCpus': 2.0,
                    'CumulativeRemoteSysCpu': 0.0,
                    'CumulativeRemoteUserCpu': 0.0,
                    'JobStatus': 2,
                    'NumShadowStarts': '1',
                    'JobRunCount': '1',
                    'RemoteHost': 'slot1_1@nid00333',
                    'QDate': '1657839749',
                },
                {
                    'ClusterId': '38612',
                    'RequestMemory': 2.0,
                    'RequestCpus': 2.0,
                    'CumulativeRemoteSysCpu': 0.0,
                    'CumulativeRemoteUserCpu': 0.0,
                    'JobStatus': 2,
                    'NumShadowStarts': '1',
                    'JobRunCount': '1',
                    'RemoteHost': 'slot1_1@nid00362',
                    'QDate': '1657839750',
                },
                {
                    'ClusterId': '39649',
                    'RequestMemory': 115.0,
                    'RequestCpus': 32.0,
                    'CumulativeRemoteSysCpu': 0.0,
                    'CumulativeRemoteUserCpu': 0.0,
                    'JobStatus': 1,
                    'NumShadowStarts': 'undefined',
                    'JobRunCount': 'undefined',
                    'RemoteHost': 'undefined',
                    'QDate': '1658136630',
                },
                {
                    'ClusterId': '39650',
                    'RequestMemory': 10.0,
                    'RequestCpus': 1.0,
                    'CumulativeRemoteSysCpu': 0.0,
                    'CumulativeRemoteUserCpu': 0.0,
                    'JobStatus': 1,
                    'NumShadowStarts': 'undefined',
                    'JobRunCount': 'undefined',
                    'RemoteHost': 'undefined',
                    'QDate': '1658145932',
                },
                {
                    'ClusterId': '39651',
                    'RequestMemory': 115.0,
                    'RequestCpus': 32.0,
                    'CumulativeRemoteSysCpu': 0.0,
                    'CumulativeRemoteUserCpu': 0.0,
                    'JobStatus': 1,
                    'NumShadowStarts': 'undefined',
                    'JobRunCount': 'undefined',
                    'RemoteHost': 'undefined',
                    'QDate': '1658152224',
                }
            ],  # noqa
         )
    ]
    for output, expected in test_data:
        result = process_condor_q(output, columns)
        # Remove total runtime since it used time.now() to calculate runtime.
        for x in result:
            x.pop('total_running_time', None)
            x.pop('cpu_percentage', None)
            x.pop('total_q_time', None)
            x.pop('JobStartDate', None)
        for x in expected:
            x.pop('JobStartDate', None)

        assert result == expected
