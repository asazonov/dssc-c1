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

        # clustering and pca
        clust_df = pd.read_csv("../data_out/clust.csv")
        pca_df = pd.read_csv("../data_out/princomps.csv")
        n_clusters = len(clust_df["cluster"].unique())

        clusters_list = list()
        pca_list = list()

        for n in range(n_clusters):
            clusters_list.append(
                {"key": "Group " + str(n + 1), "values": list()})
            pca_list.append({"key": "Group " + str(n + 1), "values": list()})

        cell_group_dict = dict()

        for i, row in clust_df.iterrows():
            cluster_cell = {
                "cell": row["cell"],
                "x": row["tsne_x"],
                "y": row["tsne_y"],
                "size": 1,
                "shape": shapes[row["cluster"] % n_clusters - 1],
                "symbol": row["cluster"],
            }
            clusters_list[row["cluster"] - 1]["values"].append(cluster_cell)
            cell_group_dict[row["cell"]] = row["cluster"]

        for i, row in pca_df.iterrows():
            pca_cell_group = cell_group_dict[row["cell"]]

            pca_cell = {
                "cell": row["cell"],
                "x": row["PC1"],
                "y": row["PC2"],
                "size": 1,
                "shape": shapes[pca_cell_group % n_clusters - 1],
                "symbol": pca_cell_group,
            }
            pca_list[pca_cell_group - 1]["values"].append(pca_cell)

        self.render("example.html",
                    vargenes_json=json.dumps(
                        [highly_variables, regular_variables]),
                    clusters_json=json.dumps(clusters_list),
                    pca_json=json.dumps(pca_list))


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
