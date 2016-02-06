from rserver import app

__all__ = ["app"]


def run(host, port):
    single_cell_rserver = app.SingleCellRServer()
    single_cell_rserver.start(host, port)
