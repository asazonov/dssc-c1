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
        fileinfo = self.request.files["file"][0]
        filename = fileinfo["filename"]
        fh = open(__UPLOADS__ + filename, 'wb+')
        fh.write(fileinfo['body'])
        self.redirect("results/" + filename.split(".")[0])


class GetResults(tornado.web.RequestHandler):
    def get(self, filename):
        self.render("done.html")
