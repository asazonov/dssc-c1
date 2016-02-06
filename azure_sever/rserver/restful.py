import tornado.web
import os
import uuid
import subprocess
import time

__UPLOADS__ = "uploads/"


# def execute(command):
#     p = subprocess.Popen(command.split(), cwd=None, stdout=subprocess.PIPE)
#     while True:
#         line = p.stdout.readline()
#         if line == '' and p.poll() is not None:
#             break
#         if line != '':
#             time_string = time.strftime("%Y-%m-%d %H:%M:%S")
#             print("\033[94m$ Rscript {time} {output}".format(time=time_string,
#                                                              output=line))
#

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
