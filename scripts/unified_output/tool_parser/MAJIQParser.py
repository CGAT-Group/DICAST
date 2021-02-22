from tool_parser.ASParser import ToolParser
from unified.Event import *
import sys, csv

#gene_id	lsv_id	num_junctions	num_exons	de_novo_junctions	seqid	strand	junctions_coords	exons_coords	ir_coords	A3SS	A5SS	ES

nan_max = sys.maxsize
nan_min = 0


class MAJIQParser(ToolParser):
    NAME = "MAJIQ"

    def __init__(self, psi_filepath, voila_filepath, outdir, combine_me):
        super().__init__("MAJIQ")
        self.combine_me = combine_me            # true -> mes, mee are counted as es events
        merged_file = merge_psi_voila(psi_filepath, voila_filepath, outdir)
        self.fh = openFile(merged_file)
        self.header = self.fh.readline().rstrip().split("\t")
        self.header_index = {}
        for i in range(len(self.header)):
            self.header_index[self.header[i]] = i
        self.ir_existent = "ir_coords" in self.header
        if self.ir_existent:
            self.header_index["ir_coords"] = self.header.index("ir_coords")
        else:
            self.header_index["ir_coords"] = -1
        self.known_event_types = ["es", "mes", "ir", "a3", "a5", "FP"]

    def closeFile(self):
        return self.fh.close()

    # return list:
    # [0]: if event is complex or not
    # [1]: what type(s) of event(s) does the event have
    def getEventType(self, event_line):
        complex_event = False
        if int(event_line[self.header_index["num_junctions"]]) > 2:  # if more than 2 junctions complex
            # TODO check this!!
            complex_event = True
        a3_bool = event_line[self.header_index["A3SS"]] == "True"
        a5_bool = event_line[self.header_index["A5SS"]] == "True"
        es_bool = event_line[self.header_index["ES"]] == "True"
        if self.ir_existent:
            ir_bool = len(event_line[self.header_index["ir_coords"]]) > 0
        else:
            ir_bool = False
        bool_sum = a3_bool + a5_bool + es_bool + ir_bool
        if bool_sum > 1:  # if more than one event type: complex
            complex_event = True
        if bool_sum == 0:
            print("No event type found for this line: " + event_line)
            return None, None
        events = []
        events.append("a3") if a3_bool else events.append("")
        events.append("a5") if a5_bool else events.append("")
        events.append("es") if es_bool else events.append("")
        events.append("ir") if ir_bool else events.append("")

        return complex_event, events

    def nextEventSet(self):
        next_line = self.fh.readline()
        if next_line == "":
            self.closeFile()
            return False
        next_line = next_line.strip('\r\n')
        events = self.parseLineToEvents(next_line)

        return events

    def parseLineToEvents(self, event_line):
        event_line = event_line.split("\t")
        gene = event_line[0]
        idx = event_line[1]
        symbol = event_line[self.header_index["seqid"]]

        strand = event_line[self.header_index["strand"]]

        junctions = event_line[self.header_index["junctions_coords"]].split(";")
        event_exons = event_line[self.header_index["exons_coords"]].split(";")
        complex_event, event_types = self.getEventType(event_line)
        # check if only IR event
        only_ir = True if not complex_event and event_types[3] == "ir" else False
        nan_fake_values = []
        if not only_ir:
            has_nan_exon, nan_exon, nan_index, nan_multiple = check_event_exons(event_exons)
            if has_nan_exon and not nan_multiple:
                event_exons, nan_fake_values = handle_nan_events(event_exons, nan_exon)
            if has_nan_exon and nan_multiple:
                event_exons, nan_fake_values = handle_nan_events_multiple(event_exons)

        if only_ir:
            ir_coords = event_line[self.header_index["ir_coords"]]
            intron_start = ir_coords.split("-")[0]
            intron_end = ir_coords.split("-")[1]
            event = IrEvent(idx, intron_start, intron_end, strand, gene, symbol)
            return [event]
        else:
            if event_types[3] == "ir":
                ir_coords = event_line[self.header_index["ir_coords"]]
                events = handle_ir_plus_complex(event_exons, junctions, event_types, idx, strand, gene, symbol, ir_coords, combine_me=self.combine_me, nan_fake_values=nan_fake_values)
            else:
                events = handle_complex_event(event_exons, junctions, event_types, idx, strand, gene, symbol, combine_me=self.combine_me, nan_fake_values=nan_fake_values)
            return events


# wrapper for all complex event types
def handle_complex_event(event_exons, junctions, event_types, idx, strand, gene, symbol, nan_fake_values, combine_me, counter_start=1):
    if event_types[0] == "a3" and event_types[1] == "a5" and event_types[2] == "":
        return handle_a3_a5_complex(event_exons, junctions, idx, strand, gene, symbol, nan_fake_values=nan_fake_values, counter_start=counter_start)
    elif event_types[0] == "a3" and event_types[1] == "a5" and event_types[2] == "es":
        out1 = set(handle_a3_a5_complex(event_exons, junctions, idx, strand, gene, symbol, nan_fake_values=nan_fake_values, counter_start=counter_start))
        out2 = set(handle_a35_es_mes_complex(event_exons, junctions, event_types, idx, strand, gene, symbol, nan_fake_values=nan_fake_values,  counter_start=counter_start, combine_me=combine_me))
        return list(out1 ^ out2)  # return symmetric difference of sets -> no duplicates
    else:
        return handle_a35_es_mes_complex(event_exons, junctions, event_types, idx, strand, gene, symbol, nan_fake_values=nan_fake_values,  counter_start=counter_start, combine_me=combine_me)


# handle a3 AND a5 events in one line
def handle_a3_a5_complex(event_exons, junctions, idx, strand, gene, symbol, nan_fake_values, counter_start=1):
    alt_events = {"a3": [], "a5": []}
    for i in range(len(junctions)):
        junc = junctions[i]
        j_start = int(junc.split("-")[0])
        j_end = int(junc.split("-")[1])
        for e in range(len(event_exons)):
            exon = event_exons[e].split("-")
            e_start, e_stop = int(exon[0]), int(exon[1])
            # exon ends before junction starts / exon starts after junction ends  -> skip exon
            if e_stop < j_start or e_start > j_end:
                continue
            # junction starts or ends inside exon -> junction is part of a3/5
            junc_start_in_exon = e_start < j_start < e_stop
            junc_end_in_exon = e_start < j_end < e_stop
            if junc_end_in_exon:
                if strand == "+":
                    alt_events["a3"].append((e_start, j_end))
                if strand == "-":
                    alt_events["a5"].append((e_start, j_end))
            if junc_start_in_exon:
                if strand == "+":
                    alt_events["a5"].append((j_start, e_stop))
                if strand == "-":
                    alt_events["a3"].append((j_start, e_stop))

    counter = counter_start
    out_events = []
    for event_type, alt_part_list in alt_events.items():
        if event_type == "a3":
            for alt_part in alt_part_list:
                a_idx = idx + ":" + str(counter)
                out_events.append(write_a35_event(event_type="a3", idx=a_idx, alt_part=alt_part, strand=strand, gene=gene, symbol=symbol,nan_fake_values=nan_fake_values, count=counter))
                counter += 1
        if event_type == "a5":
            for alt_part in alt_part_list:
                a_idx = idx + ":" + str(counter)
                out_events.append(write_a35_event(event_type="a5", idx=a_idx, alt_part=alt_part, strand=strand, gene=gene, symbol=symbol, nan_fake_values=nan_fake_values, count=counter))
                counter += 1

    return out_events


# handle IR and one other event (can be complex or basic)
def handle_ir_plus_complex(event_exons, junctions, event_types, idx, strand, gene, symbol, ir_coords, combine_me, nan_fake_values=None):
    events = []
    intron_start = ir_coords.split("-")[0]
    intron_end = ir_coords.split("-")[1]
    events.append(IrEvent(idx, intron_start, intron_end, strand, gene, symbol))
    # remove 'intron-junction' from regular junctions
    junctions.remove(ir_coords)
    remaining_events = len([0 for i in event_types if i != ''])
    if remaining_events <= 2: is_complex=False

    events.extend(handle_complex_event(event_exons, junctions, event_types, idx, strand, gene, symbol, counter_start=2, nan_fake_values=nan_fake_values, combine_me=combine_me))
    return events


# case a3/a5 & es (or mes)
def handle_a35_es_mes_complex(event_exons, junctions, event_types, idx, strand, gene, symbol, nan_fake_values, combine_me, counter_start=1):
    a3_or_a5 = "a3" if event_types[0] == "a3" else "a5"
    case1 = ((a3_or_a5 == "a3" and strand == "+") or (a3_or_a5 == "a5" and strand == "-"))
    case2 = ((a3_or_a5 == "a3" and strand == "-") or (a3_or_a5 == "a5" and strand == "+"))

    skipped_exons = {}  # store, which junction skips which exon(s)
    alt_junctions = []  # store alternative parts of exons as list of tuples
    for i in range(len(junctions)):
        junc = junctions[i]
        j_start = int(junc.split("-")[0])
        j_end = int(junc.split("-")[1])
        skipped_exons[(j_start, j_end)] = []  # init dict with junction as key -> value: skipped exon(s)
        for e in range(len(event_exons)):           # go over event exons for this junction
            if case1:
                exon = list(reversed(event_exons))[e].split("-")       # start at last exon, if looking for case1
                j_check = j_end
            elif case2:
                exon = event_exons[e].split("-")
                j_check = j_start
            else:
                return None
            e_start, e_stop = int(exon[0]), int(exon[1])
            if e_stop < j_start or e_start > j_end:         #skip exons before/after junction
                continue
            #if e_start < j_check < e_stop and not found_alt:        # exon is part of alternative a3/5; only take first exon, which fits
            if e_start < j_check < e_stop:
                if case1: alt_part = (e_start, j_end)
                if case2: alt_part = (j_start, e_stop)
                alt_junctions.append(alt_part)                      # save alt-part
                #found_alt = True
            # this exon is spanned by current junction -> part of es or mes
            if j_start < e_start and j_end > e_stop:
                skipped_exons[(j_start, j_end)].append(exon)  # save skipped exon in dict with junc as key
                found_es = True

    # generate event output objects
    out_events = []
    seen = []       # catch duplicate ES events
    counter = counter_start
    for skipped_exons_lst in skipped_exons.values():
        if len(skipped_exons_lst) == 0 or skipped_exons_lst in seen:
            continue
        seen.append(skipped_exons_lst)
        es_idx = idx + ":" + str(counter)
        mes = True if len(skipped_exons_lst) > 1 else False
        event = write_es_mes_event(idx=es_idx, strand=strand, mes=mes, skipped_exons=skipped_exons_lst, nan_fake_values=nan_fake_values, gene=gene, symbol=symbol, count=counter, combine_me=combine_me)
        out_events.extend(event)
        counter += 1
    for alt_part in alt_junctions:
        if alt_part is None:
            print("Stored alternative part of exon with None value for: " + idx)
            return None
        a_idx = idx + ":" + str(counter)
        event = write_a35_event(event_type=a3_or_a5, idx=a_idx, alt_part=alt_part, strand=strand, nan_fake_values=nan_fake_values, gene=gene, symbol=symbol, count=counter)
        out_events.append(event)
        counter += 1

    #out_events = self.rewrite_nan(out_events, nan_fake_values)

    return out_events


# function to handle basic ES events (2 junctions, 3 exons) or to simply create MES event from given skipped_exons
def write_es_mes_event(idx, strand, gene, symbol, nan_fake_values, count, combine_me, mes=False, skipped_exons=None):
    out_lst = []
    if skipped_exons is not None:
        if mes:
            mes_skipped = []
            i = 0
            for skipped in skipped_exons:
                skipped_start = "nan" if skipped[0] in nan_fake_values else str(skipped[0])
                skipped_end = "nan" if skipped[1] in nan_fake_values else str(skipped[1])
                mes_skipped.append([skipped_start, skipped_end])
                if combine_me:
                    out_lst.append(EsEvent(idx, skipped_start, skipped_end, strand, gene, symbol=symbol, count=count+i))
                i += 1
            if not combine_me:
                out_lst.append(MesEvent(idx, mes_skipped=mes_skipped, strand=strand, gene=gene, symbol=symbol, count=count))
        else:
            skipped_start = "nan" if skipped_exons[0][0] in nan_fake_values else skipped_exons[0][0]
            skipped_end = "nan" if skipped_exons[0][1] in nan_fake_values else skipped_exons[0][1]
            out_lst.append(EsEvent(idx, skipped_start, skipped_end, strand, gene, symbol=symbol, count=count))
        return out_lst
    else:
        print("No skipped exons given:" + idx)
        return [None]


# create A3/5 event by giving the alternative region in alt_part
def write_a35_event(event_type, idx, strand, gene, symbol, nan_fake_values, count, alt_part=None):
    if alt_part is not None:
        alt1 = "nan" if alt_part[0] in nan_fake_values else alt_part[0]
        alt2 = "nan" if alt_part[1] in nan_fake_values else alt_part[1]
        if event_type == "a3":
            return A3Event(idx, alt1, alt2, strand, gene, symbol, count=count)
        if event_type == "a5":
            return A5Event(idx, alt1, alt2, strand, gene, symbol, count=count)
    else:
        print("No alternative part given:" + idx)
        return None


# function to check the event_exons:
# if NaN in one of coords -> exon wont be in gtf file -> regular functions will not work
# check if there are gtf_exons in between event_exons (to check for MES event)
def check_event_exons(event_exons):
    has_nan_exon = False
    nan_exon = None
    nan_index = None
    nan_count = 0

    for i in range(len(event_exons)):
        ee = event_exons[i]
        e_split = ee.split("-")
        ee_start = e_split[0]
        ee_end = e_split[1]
        if ee_start == "nan" or ee_end == "nan":
            has_nan_exon = True
            nan_exon = ee
            nan_index = i
            nan_count += 1

    nan_multiple = nan_count > 1
    return has_nan_exon, nan_exon, nan_index, nan_multiple


def handle_nan_events_multiple(event_exons):
    event_exons_copy = []
    nan_fake_values = []
    for e in range(len(event_exons)):
        exon = event_exons[e].split("-")
        if exon[0] != "nan" and exon[1] != "nan":
            event_exons_copy.append(event_exons[e])
            continue
        elif exon[0] == "nan":
            coord = int(exon[1])-1
            event_exons_copy.append(str(coord)+"-"+exon[1])
            nan_fake_values.append(str(coord))
        elif exon[1] == "nan":
            coord = int(exon[0])+1
            event_exons_copy.append(exon[0]+"-"+str(coord))
            nan_fake_values.append(str(coord))
    return event_exons_copy, nan_fake_values


def handle_nan_events(event_exons, nan_exon):
    event_exons_copy = event_exons.copy()
    event_exons_copy.remove(nan_exon)
    nan_exon = nan_exon.split("-")
    if nan_exon[0] == "nan":
        nan_coordinate = int(nan_exon[1])
        nan_fake_value = str(nan_coordinate-1)
        insert_exon = str(nan_fake_value)+"-"+str(nan_coordinate)
    else:
        nan_coordinate = int(nan_exon[0])
        nan_fake_value = str(nan_coordinate+1)
        insert_exon = str(nan_coordinate)+"-"+str(nan_fake_value)

    for i in range(len(event_exons_copy) - 1):
        exon = event_exons_copy[i].split("-")
        e_stop = int(exon[1])
        exon_next = event_exons_copy[i + 1].split("-")
        e_next_start = int(exon_next[0])
        if e_stop < nan_coordinate < e_next_start:
            event_exons_copy.insert(i+1, insert_exon)
            break

    return event_exons_copy, [nan_fake_value]


def merge_psi_voila(psi_file, voila_file, outdir):
    merged_file = outdir+"/merged_majiq.tsv"
    header = ["gene_id", "lsv_id", "num_junctions", "num_exons", "de_novo_junctions", "seqid", "strand",
              "junctions_coords", "exons_coords", "ir_coords", "A3SS", "A5SS", "ES"]

    psi_dict = {}
    with open(psi_file, "r") as f:
        # read header
        f.readline()
        for line in f:
            line = line.strip().split("\t")
            idx = line[1]
            a5 = line[5]
            a3 = line[6]
            es = line[7]
            psi_dict[idx] = [a3, a5, es]
    f.close()

    voila_dict = {}
    with open(voila_file, 'r') as f:
        while True:
            line = f.readline()
            if not line.startswith("#"): break
        h = line.split("\t")
        for line in f:
            line = line.strip().split("\t")
            idx = line[1]
            voila_dict[idx] = [line[0], line[1], line[5], line[6], line[7], line[8], line[9], line[10], line[11],
                               line[12]]
    # gene_id	lsv_id	mean_psi_per_lsv_junction	stdev_psi_per_lsv_junction	lsv_type	num_junctions	num_exons	de_novo_junctions	seqid	strand	junctions_coords	exons_coords	ir_coords	ucsc_lsv_link
    # 0           1                                                                                   5             6             7                8       9          10               11            12
    with open(merged_file, 'w') as out:
        tsv_writer = csv.writer(out, delimiter='\t')
        tsv_writer.writerow(header)
        for idx, values in voila_dict.items():
            values.extend(psi_dict[idx])
            tsv_writer.writerow(values)

    return merged_file


def openFile(filepath):
    return open(filepath, mode='r')