import tornado.web
import os
import uuid
import subprocess

__UPLOADS__ = "uploads/"


class Index(tornado.web.RequestHandler):
    """
    Returns the index page.
    """

    def get(self):
        subprocess.call(["Rscript", "../helloworld.R"])
        self.render("index.html")


class UploadData(tornado.web.RequestHandler):
    def post(self):
        # save file
        fileinfo = self.request.files['filearg'][0]
        print("fileinfo is", fileinfo)
        fname = fileinfo['filename']
        extn = os.path.splitext(fname)[1]
        cname = str(uuid.uuid4()) + extn
        fh = open(__UPLOADS__ + cname, 'w')
        fh.write(fileinfo['body'])
        self.render("done.html")


class GetResults(tornado.web.RequestHandler):
    def get(self):
        try:
            result = self.get_argument("result")
            self.application.broadcast_tracking_result(result)
            self.finish({"response": "OK"})
        except tornado.web.MissingArgumentError:
            self.set_status(400)
            self.finish("???")
