import re, argparse, dataclasses

from .build     import get_mfc_target_names, get_target_names, get_dependencies_names
from .common    import format_list_to_string
from .test.test import CASES as TEST_CASES


def parse(config):
    from .run.engines  import ENGINES
    from .run.mpi_bins import BINARIES

    parser = argparse.ArgumentParser(
        prog="./mfc.sh",
        description=f"""\
Welcome to the MFC master script. This tool automates and manages building, testing, \
running, and cleaning of MFC in various configurations on all supported platforms. \
The README documents this tool and its various commands in more detail. To get \
started, run ./mfc.sh build -h.""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parsers = parser.add_subparsers(dest="command")

    run   = parsers.add_parser(name="run",   help="Run a case with MFC.",            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    test  = parsers.add_parser(name="test",  help="Run MFC's test suite.",           formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    build = parsers.add_parser(name="build", help="Build MFC and its dependencies.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    clean = parsers.add_parser(name="clean", help="Clean build artifacts.",          formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    bench = parsers.add_parser(name="bench", help="Benchmark MFC (for CI).",         formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    count = parsers.add_parser(name="count", help="Count LOC in MFC.",         formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    def add_common_arguments(p, mask = None):
        if mask is None:
            mask = ""

        if "t" not in mask:
            p.add_argument("-t", "--targets", metavar="TARGET", nargs="+", type=str.lower, choices=get_target_names(),
                           default=get_mfc_target_names(), help=f"Space separated list of targets to act upon. Allowed values are: {format_list_to_string(get_target_names())}.")

        if "m" not in mask:
            for f in dataclasses.fields(config):
                p.add_argument(   f"--{f.name}", action="store_true",               help=f"Turn the {f.name} option ON.")
                p.add_argument(f"--no-{f.name}", action="store_false", dest=f.name, help=f"Turn the {f.name} option OFF.")

            p.set_defaults(**{ f.name: getattr(config, f.name) for f in dataclasses.fields(config) })
            
        if "j" not in mask:
            p.add_argument("-j", "--jobs", metavar="JOBS", type=int, default=1, help="Allows for JOBS concurrent jobs.")

        if "v" not in mask:
            p.add_argument("-v", "--verbose", action="store_true", help="Enables verbose compiler & linker output.")

        if "n" not in mask:
            for name in get_dependencies_names():
                p.add_argument(f"--no-{name}", action="store_true", help=f"Do not build the {name} dependency. Use the system's instead.")

    # === BUILD ===
    add_common_arguments(build)

    # === CLEAN ===
    add_common_arguments(clean, "j")

    binaries = [ b.bin for b in BINARIES ]

    # === TEST ===
    add_common_arguments(test, "t")
    test.add_argument(      "--generate",    action="store_true", help="Generate golden files.")
    test.add_argument("-l", "--list",        action="store_true", help="List all available tests.")
    test.add_argument("-f", "--from",        default=TEST_CASES[0].get_uuid(), type=str, help="First test UUID to run.")
    test.add_argument("-t", "--to",          default=TEST_CASES[-1].get_uuid(), type=str, help="Last test UUID to run.")
    test.add_argument("-o", "--only",        nargs="+", type=str, default=[], metavar="L", help="Only run tests with UUIDs or hashes L.")
    test.add_argument("-b", "--binary",      choices=binaries, type=str, default=None, help="(Serial) Override MPI execution binary")
    test.add_argument("-r", "--relentless",  action="store_true", default=False, help="Run all tests, even if multiple fail.")
    test.add_argument("-a", "--test-all",    action="store_true", default=False, help="Run the Post Process Tests too.")
    test.add_argument("--case-optimization", action="store_true", default=False, help="(GPU Optimization) Compile MFC targets with some case parameters hard-coded.")

    # === RUN ===
    engines = [ e.slug for e in ENGINES ]

    add_common_arguments(run)
    run.add_argument("input",                  metavar="INPUT",                 type=str,                     help="Input file to run.")
    run.add_argument("arguments",              metavar="ARGUMENTS", nargs='*',  type=str, default=[],         help="Additional arguments to pass to the case file.")
    run.add_argument("-e", "--engine",         choices=engines,                 type=str, default=engines[0], help="Job execution/submission engine choice.")
    run.add_argument("-p", "--partition",      metavar="PARTITION",             type=str, default="",         help="(Batch) Partition for job submission.")
    run.add_argument("-N", "--nodes",          metavar="NODES",                 type=int, default=1,          help="(Batch) Number of nodes.")
    run.add_argument("-n", "--tasks-per-node", metavar="TASKS",                 type=int, default=1,          help="Number of tasks per node.")
    run.add_argument("-w", "--walltime",       metavar="WALLTIME",              type=str, default="01:00:00", help="(Batch) Walltime.")
    run.add_argument("-a", "--account",        metavar="ACCOUNT",               type=str, default="",         help="(Batch) Account to charge.")
    run.add_argument("-@", "--email",          metavar="EMAIL",                 type=str, default="",         help="(Batch) Email for job notification.")
    run.add_argument("-#", "--name",           metavar="NAME",                  type=str, default="MFC",      help="(Batch) Job name.")
    run.add_argument("-f", "--flags",          metavar="FLAGS",     nargs='+',  type=str, default=[],         help="(Batch) Additional batch options.")
    run.add_argument("-b", "--binary",         choices=binaries,                type=str, default=None,       help="(Interactive) Override MPI execution binary")
    run.add_argument("-s", "--scratch",        action="store_true",                       default=False,      help="Build from scratch.")
    run.add_argument("--ncu",                  nargs=argparse.REMAINDER,        type=str,                     help="Profile with NVIDIA Nsight Compute.")
    run.add_argument("--nsys",                 nargs=argparse.REMAINDER,        type=str,                     help="Profile with NVIDIA Nsight Systems.")
    run.add_argument(      "--dry-run",        action="store_true",                       default=False,      help="(Batch) Run without submitting batch file.")
    run.add_argument("--case-optimization",    action="store_true",                       default=False,      help="(GPU Optimization) Compile MFC targets with some case parameters hard-coded.")
    run.add_argument(      "--no-build",       action="store_true",                       default=False,      help="(Testing) Do not rebuild MFC.")
    run.add_argument("--wait",                 action="store_true",                       default=False,      help="(Batch) Wait for the job to finish.")

    # === BENCH ===
    add_common_arguments(bench, "t")

    # === COUNT ===
    add_common_arguments(count)

    args: dict = vars(parser.parse_args())

    # Add default arguments of other subparsers
    def append_defaults_to_data(name: str, parser):
        if args["command"] != name:
            vals, errs = parser.parse_known_args(["-i None"])
            for key,val in vars(vals).items():
                if not key in args:
                    args[key] = val

    for a, b in [("run",   run  ), ("test",  test ), ("build", build),
                 ("clean", clean), ("bench", bench), ("count", count)]:
        append_defaults_to_data(a, b)

    if args["command"] is None:
        parser.print_help()
        exit(-1)

    # "Slugify" the name of the job
    args["name"] = re.sub(r'[\W_]+', '-', args["name"])

    return args
