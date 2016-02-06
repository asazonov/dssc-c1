import logging
import os
import signal
import tornado.httpserver
import tornado.ioloop
import tornado.log
import tornado.web
from rserver import urls
from rserver import restful


class SingleCellRServer(tornado.web.Application):
    def __init__(self):
        handlers = [
            tornado.web.url(urls.index, restful.Index),
            tornado.web.url(urls.post_upload, restful.UploadData),
            tornado.web.url(urls.get_result, restful.GetResults),
        ]

        settings = {
            "static_path": os.path.join(os.path.dirname(__file__), "static"),
            "template_path": os.path.join(os.path.dirname(__file__),
                                          "templates"),
            "autoescape": None,
            "autoreload": False,
            "debug": True,
        }

        tornado.web.Application.__init__(self, handlers, **settings)

        tornado.log.enable_pretty_logging()

        self.http_server = tornado.httpserver.HTTPServer(self)

        self.is_closing = False

    def start(self, host, port):
        self.http_server.bind(address=host, port=port)
        self.http_server.start()
        logging.info("Starting @ {addr}:{port}".format(addr=host, port=port))
        signal.signal(signal.SIGINT, self.signal_handler)
        tornado.ioloop.PeriodicCallback(self.try_exit, 100).start()
        tornado.ioloop.IOLoop.instance().start()

    def signal_handler(self, signum, frame):
        self.is_closing = True
        logging.info("Exiting...")

    def try_exit(self):
        if self.is_closing:
            tornado.ioloop.IOLoop.instance().stop()
            logging.info("Exit success")
