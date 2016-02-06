import rserver
import argparse


def main():
    parser = argparse.ArgumentParser(
        description="Run the Azure R server"
    )

    parser.add_argument(
        "--host", dest="host", type=str, default="0.0.0.0",
        help="IP address of the R server"
    )

    parser.add_argument(
        "--port", dest="port", type=int, default=5000,
        help="Port of the R server"
    )

    args = parser.parse_args()

    rserver.run(args.host, args.port)


if __name__ == "__main__":
    main()
