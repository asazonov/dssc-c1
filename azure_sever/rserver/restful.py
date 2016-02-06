import tornado.web
import os
import uuid

__UPLOADS__ = "uploads/"


class Index(tornado.web.RequestHandler):
    """
    Returns the index page.
    """

    def get(self):
        self.render("index.html")


class UploadData(tornado.web.RequestHandler):
    def post(self):
        fileinfo = self.request.files['filearg'][0]
        fname = fileinfo['filename']
        extn = os.path.splitext(fname)[1]
        cname = str(uuid.uuid4()) + extn
        fh = open(__UPLOADS__ + cname, 'w')
        fh.write(fileinfo['body'])
        self.finish(cname + " is uploaded!! Check %s folder" % __UPLOADS__)


class GetResults(tornado.web.RequestHandler):
    def get(self):
        try:
            result = self.get_argument("result")
            self.application.broadcast_tracking_result(result)
            self.finish({"response": "OK"})
        except tornado.web.MissingArgumentError:
            self.set_status(400)
            self.finish("???")
