from jaws_jtm.lib.jtminterface import JtmInterface


def check_manager():
    desc = u"jtm-check-manager"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-cl", "--cluster",
                        dest="cluster_name",
                        required=False,
                        choices=COMPUTE_RESOURCES)
    parser.add_argument("-l", "--loglevel",
                        help="Set loglevel (default=info).",
                        dest="log_level", required=False,
                        default=None)
    args = parser.parse_args()

    add_info = socket.gethostname().replace(".", "_")
    if args.cluster_name:
        add_info = args.cluster_name

    ret = int(JtmInterface('check_manager',
                           info_tag=add_info)).call(jtm_host_name=args.cluster_name, log_level=args.log_level))
    if ret is None or ret == -88:
        sys.exit(-1)
    sys.exit(0) if ret != 0 else sys.exit(1)
