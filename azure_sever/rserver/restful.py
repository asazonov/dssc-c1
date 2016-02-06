import tornado.web


class Index(tornado.web.RequestHandler):
    """
    Returns the index page.
    """

    def get(self):
        self.render("index.html")


class UploadData(tornado.web.RequestHandler):
    def post(self):
        try:
            result = self.get_argument("result")
            self.application.broadcast_tracking_result(result)
            self.finish({"response": "OK"})
        except tornado.web.MissingArgumentError:
            self.set_status(400)
            self.finish("???")


class GetResults(tornado.web.RequestHandler):
    def get(self):
        try:
            result = self.get_argument("result")
            self.application.broadcast_tracking_result(result)
            self.finish({"response": "OK"})
        except tornado.web.MissingArgumentError:
            self.set_status(400)
            self.finish("???")
