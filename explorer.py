import os
import numpy
import re
import scipy
import scipy.stats
import math
import plotly.offline as off
off.init_notebook_mode(connected=True)
import plotly.graph_objs as go
import matplotlib.pyplot as plt

# copied from predictions.ipynb
class Metadata:
    def __init__(self, ids, cels, metadata):
        self.IDs = ids
        self.CELs = cels
        self.modal_allele = [int(metadata[i]["modal_allele"]) for i in ids]
        self.progenitor_allele = [int(metadata[i]["progenitor_allele"]) for i in ids]
        self.MIRS = [int(metadata[i]["MIRS"]) for i in ids]
    def __str__(self):
        return "Metadata={{IDs: {}...,\n CELs: {}...,\n modal_allele: {}...,\n progenitor_allele: {}...,\n MIRS: {}...}}".format(self.IDs[:5], self.CELs[:5], self.modal_allele[:5], self.progenitor_allele[:5], self.MIRS[:5])
    def __repr__(self):
        return self.__str__()

def load_metadata(metadata_path):
    metadata = {}
    metadata_order = []
    with open(metadata_path) as f:
        for i, line in enumerate(f):
            line = line.strip().split()
            if i == 0:
                names = line[1:]
            else:
                values = line[1:]
                patient_id = line[0]
                metadata_order.append(patient_id)
                metadata[patient_id] = {k: v for k, v in zip(names, values)}
    blood_IDs = [i for i in metadata_order]
    muscle_IDs = [i for i in metadata_order if metadata[i]["muscle_cel"] != "refused_biopsy"]
    blood_CELs = [metadata[i]["blood_cel"] for i in blood_IDs]
    muscle_CELs = [metadata[i]["muscle_cel"] for i in muscle_IDs]
        
    blood_record = Metadata(blood_IDs, blood_CELs, metadata)
    muscle_record = Metadata(muscle_IDs, muscle_CELs, metadata)
    return blood_record, muscle_record

# copied from predictions.ipynb
def load_genes(genecode_genes, filename):
    genes = set()
    repeated = set()
    with open(filename) as f:
        for line in f:
            line = line.rstrip()
            if line in genes:
                repeated.add(line)
            genes.add(line)
    failed = genes.difference(genecode_genes)
    if repeated:
        print("Loading {}. These genes appear more than once: {}".format(filename, [g for g in repeated]))
    if failed:
        print("Loading {}. Couldn't identify the following genes: {}".format(filename, [g for g in failed]))
    return genes.intersection(genecode_genes)

# copied from predictions.ipynb
def produce_data(meta, gene_names, experiment):
    IDs = meta.IDs
    probe_data = []
    probe_IDs = []
    for gene in gene_names:        
        with open(os.path.join(experiment, gene)) as f:
            for i, line in enumerate(f):
                line = line.rstrip().split()
                if i == 0:
                    prefix = "patient_"
                    our_IDs = [elem[len(prefix):] for elem in line if re.match(prefix, elem)]
                    try:
                        assert IDs == our_IDs
                    except AssertionError:
                        print(gene)
                    headers = {header: i for i, header in enumerate(line)}
                    patient_data = {header[len(prefix):]: i for i, header in enumerate(line) if re.match(prefix, header)}
                    def write_signature(line):
                        signature = []
                        for elem in ["gene_name", "probeset_id", "seq5to3plus", "chrom", "strand", "genocode_left", "genecode_right"]:
                            signature.append(line[headers[elem]])
                        return signature
                else:
                    probe_ID = write_signature(line)
                    rv = []
                    for patient_id in IDs:
                        rv.append(float(line[patient_data[patient_id]]))
                    probe_data.append(rv)
                    probe_IDs.append(probe_ID)
    probe_data = numpy.array(probe_data)
    return probe_data, probe_IDs

def name_to_plot(probesets, gene_name, regressor):
    
    slope_pvalue_intensity_position = []
    for probeset_id, probeset_data in probesets.items():
        left = probeset_data[0]
        probeset_expression = probeset_data[1]
        for probeset in probeset_expression:
            slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(regressor, probeset)
            slope_pvalue_intensity_position.append([slope, p_value, sum(probeset)/len(probeset), left, probeset_id])

    zipped_properties = list(zip(*slope_pvalue_intensity_position))
    slope, intensity, pvalues, position, probeset_belongs = zipped_properties[0], zipped_properties[2], zipped_properties[1], zipped_properties[3], zipped_properties[4]
    signed_pvalues = [-math.log(pvalue, 10) * (s/abs(s)) for pvalue, s in zip(pvalues, slope)]
    
    data = [
    go.Scatter(
      x = position,
      y = [-math.log(0.05, 10)] * len(position),
      marker = dict(
        color = 'rgba(255, 0, 0, 0.4)',
        line = dict(
            color = 'rgba(255, 0, 0, 0.4)'
        )
      ),
      mode = 'lines',
      name = 'spliced-in @ 0.05',
    ),
    go.Scatter(
      x = position,
      y = [-math.log(0.05/len(position), 10)] * len(position),
      marker = dict(
        color = 'rgba(255, 0, 0, 1)',
        line = dict(
            color = 'rgba(255, 0, 0, 1)'
        )
      ),
      mode = 'lines',
      name = 'spliced-in @ Bonferroni 0.05'
    ),
    go.Scatter(
      x = position,
      y = [math.log(0.05, 10)] * len(position),
      marker = dict(
        color = 'rgba(255, 0, 0, 0.4)',
        line = dict(
            color = 'rgba(255, 0, 0, 0.4)'
        )
      ),
      mode = 'lines',
      name = 'spliced-out @ 0.05'
    ),
    go.Scatter(
      x = position,
      y = [math.log(0.05/len(position), 10)] * len(position),
      marker = dict(
        color = 'rgba(255, 0, 0, 1)',
        line = dict(
            color = 'rgba(255, 0, 0, 1)'
        )
      ),
      mode = 'lines',
      name = 'spliced-out @ Bonferroni 0.05'
    ),
    dict(
      type = 'scatter',
      x = position,
      y = signed_pvalues,
      mode = 'markers',
      transforms = [dict(
        type = 'groupby',
        groups = probeset_belongs
      )]
    )
    ]

    #off.iplot({'data': data}, validate=False)
    off.iplot({'data': data, 'layout': dict(title = gene_name)}, validate=False)
    
def plot_gene(gene_name, metadata, location, regressor):
    expression, probe_data = produce_data(metadata, [gene_name], location)
    probesets = {}
    for i, (_, probeset_id, _, _, _, left, _) in enumerate(probe_data):
        expression_values = probesets.setdefault(probeset_id, [left, []])
        expression_values[1].append(expression[i,:])
    name_to_plot(probesets, gene_name, regressor)
    
def expression_per_probeset(probeset_id, genes, experiment, muscle_metadata):
    expression, probe_data = produce_data(muscle_metadata, genes, experiment)

    expression_indices = [i for i, elem in enumerate(probe_data) if elem[1] == probeset_id]

    data_to_plot = expression[expression_indices,:]

    for i in range(len(expression_indices)):
        plt.figure()
        x = muscle_metadata.modal_allele
        y = data_to_plot[i]
        plt.scatter(x, y, marker=",")
        plt.plot(numpy.unique(x), numpy.poly1d(numpy.polyfit(x, y, 1))(numpy.unique(x)), linewidth=4, color="r")
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
        print("p_value, slope, r_value, std_err", p_value, slope, r_value, std_err)
        plt.xlabel("Measured allele length")
        plt.ylabel("Expression level of a probe belonging to probeset with id {}".format(probeset_id))