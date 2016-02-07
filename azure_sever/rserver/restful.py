import tornado.web
import pandas as pd
import json
import math
import random
import os
import uuid
import subprocess

__UPLOADS__ = "uploads/"

shapes = ["circle", "cross", "triangle-up", "triangle-down", "diamond",
          "square"]

colors = []


class Index(tornado.web.RequestHandler):
    """
    Returns the index page.
    """

    def get(self):
        # subprocess.call(["Rscript", "../helloworld.R"])

        # significance
        highly_variables = {"key": "Highly Variable", "values": list()}
        regular_variables = {"key": "Regular", "values": list()}

        vargenes_df = pd.read_csv("../data_out/vargenes.csv")
        for i, row in vargenes_df.iterrows():
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

        # take random 800 samples
        highly_variables["values"] = random.sample(highly_variables["values"],
                                                   800)
        regular_variables["values"] = random.sample(regular_variables["values"],
                                                    800)

        # clustering
        # highly_variables = {"key": "Highly Variable", "values": list()}
        # regular_variables = {"key": "Regular", "values": list()}

        clust_df = pd.read_csv("../data_out/clust.csv")
        n_clusters = len(clust_df["cluster"].unique())

        clusters_list = list()
        for n in range(n_clusters):
            clusters_list.append(
                {"key": "Group " + str(n + 1), "values": list()})

        for i, row in clust_df.iterrows():
            gene = {
                "cell": row["cell"],
                "x": row["tsne_x"],
                "y": row["tsne_y"],
                "size": 1,
                "shape": shapes[row["cluster"] % n_clusters - 1],
                "symbol": row["cluster"],
            }

            clusters_list[row["cluster"] - 1]["values"].append(gene)

        self.render("example.html",
                    vargenes_json=json.dumps(
                        [highly_variables, regular_variables]),
                    clusters_json=json.dumps(clusters_list))


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
