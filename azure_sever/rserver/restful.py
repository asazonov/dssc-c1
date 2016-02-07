import tornado.web
import pandas as pd
import json
import math
import random
import os
import uuid
import subprocess

__UPLOADS__ = "uploads/"

shapes = ["thin-x", "circle", "cross", "triangle-up", "triangle-down",
          "diamond", "square"]

colors = []


class Index(tornado.web.RequestHandler):
    """
    Returns the index page.
    """

    def get(self):
        # subprocess.call(["Rscript", "../helloworld.R"])

        highly_variables = {"key": "Highly Variable", "values": list()}
        regular_variables = {"key": "Regular", "values": list()}

        vargenes_df = pd.read_csv("../data_out/vargenes.csv")
        for i, row in vargenes_df.iterrows():
            if pd.isnull(row["counts.avg"]) or pd.isnull(
                    row["counts.cv2"]) or pd.isnull(row["signif_var"]):
                continue

            gene = {
                "gene": row["gene"],
                "x": math.log(row["counts.avg"]),
                "y": math.log(row["counts.cv2"]),
                "size": 1
            }

            if row["signif_var"] == 1:
                gene["shape"] = "triangle-up"
                gene["symbol"] = 1
                highly_variables["values"].append(gene)
            else:
                gene["shape"] = "cross"
                gene["symbol"] = 0
                regular_variables["values"].append(gene)

        # take random 500 samples
        highly_variables["values"] = random.sample(highly_variables["values"],
                                                   800)
        regular_variables["values"] = random.sample(regular_variables["values"],
                                                    800)

        self.render("example.html", vargenes_json=json.dumps(
            [highly_variables, regular_variables]))


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
