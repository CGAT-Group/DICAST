from unified.EventAnnotationReader import EventAnnotationReader
from unified.Event import *
import csv


#chr	gene	id	strand	event_type	count	start_coordinates	end_coordinates

EVENT_TYPES = {"ES": EsEvent, "IR": IrEvent, "A3": A3Event, "A5": A5Event,
               "MES": MesEvent, "MEE": MeeEvent, "ALE": AleEvent, "AFE": AfeEvent}


def compare_with_annotation(anno: EventAnnotationReader, file, outfile, strict, threshold):
    events = anno.events
    found_events = {"ES": 0, "IR": 0, "A3": 0, "A5": 0, "MES": 0, "MEE": 0, "ALE": 0, "AFE": 0}
    correct_events = {"ES": 0, "IR": 0, "A3": 0, "A5": 0, "MES": 0, "MEE": 0, "ALE": 0, "AFE": 0}
    anno_events = anno.count_event_types
    cutoff = int(threshold)
    in_cutuff = {"ES": 0, "IR": 0, "A3": 0, "A5": 0, "MES": 0, "MEE": 0, "ALE": 0, "AFE": 0}
    seen_es = set()

    with open(file, 'r') as f:
        f.readline() # read header

        for line in f:
            event = line.strip("\n").split("\t")
            gene = event[1]
            strand = event[3]
            count = event[5]
            event_type = event[4]
            found_events[event_type] += 1   # this line could be pushed under the if check..
            # get all events from this gene with same event_type
            if gene not in events.keys():
                continue
            # get all events from gene in anno with same event type as event_to_check
            anno_events_same_type = [ev for ev in events[gene] if type(ev) is EVENT_TYPES[event_type]]
            #if event_type == "ES" and compare_es_with_me:
            #    anno_events_same_type = [ev for ev in events[gene] if type(ev) in [EVENT_TYPES["ES"], EVENT_TYPES["MES"], EVENT_TYPES["MEE"]]]
            for anno_event in anno_events_same_type:
                if event_type != "MES" and event_type != "MEE":
                    st_st = anno_event.get_start_stop()
                    if event_type == "ES":
                        seen_es.add(event[2])

                    start_anno = st_st[0]
                    stop_anno = st_st[1]

                    start_event = event[6]
                    stop_event = event[7]

                    if equal(start_event, stop_event, start_anno, stop_anno, pseudo=False, strict=strict):
                        correct_events[event_type] += 1
                        break
                    else:
                        if total_distance(start_event, stop_event, start_anno, stop_anno) <= cutoff:
                            correct_events[event_type] += 1
                            in_cutuff[event_type] += 1
                            break
                else:
                    anno_exons = anno_event.get_start_stop()
                    event_starts = event[6].split(",")
                    event_stops = event[7].split(",")

                    anno_starts = [i[0] for i in anno_exons]
                    anno_stops = [i[1] for i in anno_exons]

                    if equal_multiple(event_starts, event_stops, anno_starts, anno_stops, pseudo=False, strict=strict):
                        correct_events[event_type] += 1
                        break

    f.close()

    precision = scores_per_type(correct_events, found_events)
    recall = scores_per_type(correct_events, anno_events)

    if outfile is None:
        print("seen es:"+str(len(seen_es)))
        print("#######################################")
        print("correct_events: " + str(correct_events))
        print("found_events: " + str(found_events))
        print("anno events: " + str(anno_events))
        print("in cutoff: " + str(in_cutuff))
        print("#######################################")

        print("Precision: " + str(precision))
        print("----------------")
        print("Recall: " + str(recall))
    else:
        with open(outfile, "w") as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(["type", "correct_events", "found_events", "annotation_events", "precision", "recall"])

            for ev_type in EVENT_TYPES.keys():
                writer.writerow([ev_type, correct_events[ev_type], found_events[ev_type], anno_events[ev_type], precision[ev_type], recall[ev_type]])

        f.close()


def total_distance(start1, stop1, start2, stop2):
    if start1 == "nan" or start2 == "nan" or stop1 == "nan" or stop2 == "nan":
        diff = 999999
    else:
        diff = abs(int(start1) - int(start2)) + abs(int(stop1) -  int(stop2))
    return diff


# calculate score per event type
def scores_per_type(correct, all_found):
    result = {}
    for key in correct:
        if key in all_found:
            if all_found[key] > 0:
                result[key] = round(correct[key] / all_found[key], 3)
            else:
                result[key] = 0
        else:
            pass
    return result


def reverse_equal(start_event, stop_event, start_anno, stop_anno):
    return equal(stop_event, start_event, start_anno, stop_anno)


# compare event with event in annotation (only use on events with single start/stop)
# use with: ES, IR, A3, A5, ALE, AFE
# pseudo: True if allow for shift of +/-1 of event when comparing to anno
# strict: True if both start & stop have to be equal
def equal(start_event, stop_event, start_anno, stop_anno, pseudo=False, strict=False):

    coords = handle_nan_in_list([start_event, stop_event, start_anno, stop_anno])
    start_event = int(coords[0])
    stop_event = int(coords[1])
    start_anno = int(coords[2])
    stop_anno = int(coords[3])

    if pseudo and not strict:
        if start_event+1 == start_anno or stop_event+1 == stop_anno:
            return True
        if start_event-1 == start_anno or stop_event-1 == stop_anno:
            return True
    if pseudo and strict:
        if start_event+1 == start_anno and stop_event+1 == stop_anno:
            return True
        if start_event-1 == start_anno and stop_event-1 == stop_anno:
            return True
    if not pseudo and strict:
        if start_event == start_anno and stop_event == stop_anno:
            return True
    if not pseudo and not strict:
        if start_event == start_anno or stop_event == stop_anno:
            return True

    return False


# compare MEE and MES events with annotation
def equal_multiple(event_starts, event_stops, anno_starts, anno_stops, pseudo=False, strict=False):
    # not equal number of skipped exons -> not correct
    if len(anno_starts) != len(event_starts) or len(anno_stops) != len(event_stops):
        return False

    # check each list for nan values -> set to 0 for comparison
    coordinate_lists = [handle_nan_in_list(i) for i in [event_starts, event_stops, anno_starts, anno_stops]]

    # check for each start/stop, if it exists in the annotation

    for i in range(len(event_starts)):
        start_event = int(event_starts[i])
        stop_event = int(event_stops[i])
        correct = False
        for j in range(len(anno_starts)):
            start_anno = int(anno_starts[j])
            stop_anno = int(anno_stops[j])

            if pseudo and not strict:
                if start_event + 1 == start_anno or stop_event + 1 == stop_anno:
                    correct = True
                if start_event - 1 == start_anno or stop_event - 1 == stop_anno:
                    correct = True
            if pseudo and strict:
                if start_event + 1 == start_anno and stop_event + 1 == stop_anno:
                    correct = True
                if start_event - 1 == start_anno and stop_event - 1 == stop_anno:
                    correct = True
            if not pseudo and strict:
                if start_event == start_anno and stop_event == stop_anno:
                    correct = True
            if not pseudo and not strict:
                if start_event == start_anno or stop_event == stop_anno:
                    correct = True

        # did not find any match for this exon in the annotation -> whole MEE/MES is not correct
        if not correct:
            return False

    return True


# check for nan-values -> set nan to 0 for comparison
def handle_nan_in_list(l):
    for index, val in enumerate(l):
        if val == 'nan':
            l[index] = 0
    return l