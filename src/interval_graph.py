import matplotlib.pyplot as plt


file_path = "r4.fastq"
intervals = [] # intervals contains interval of the form [startpos, endpos, length, id]


def is_not_maximal(l, s, id):
    for [sp, ep, _, idp] in intervals:
        if idp != id and sp <= s and s+l <= ep:
            return True
    return False

def clean():
    cleaned_intervals = []
    for [s,e,l,id] in intervals:
        if is_not_maximal(l,s, id) == False:
            cleaned_intervals.append([s,e,l,id])
    return cleaned_intervals
            

# insert to make intervals sorted by increasing startpos
def insert(l, s, id):
    for i in range(len(intervals)):
        if intervals[i][0] > s:
            intervals.insert(i, [s,s+l, l, id])
            return
    intervals.append([s,s+l, l, id])


with open(file_path, "r") as f:
    for line in f.readlines():
        if line.startswith("@"):
            for s in line.split(";"):
                if s.startswith("length="):
                    length = int(s[7:-2])
                if s.startswith("startpos="):
                    startpos = int(s[9:])
                if s.startswith("@m"):
                    for s2 in s.split("/"):
                        if s2.startswith("@m"):
                            id = int(s2[2:])
            print(f"{id} {startpos} {startpos+length} {length}")
            insert(length, startpos, id)
        # if line[0] == "A" or line[0] == "C" or line[0] == "G" or line[0] == "G":
        #     print(line[:10])

    intervals = clean()    

    print(intervals)
    fig, ax = plt.subplots(figsize=(10, 6))

    rank = []
    i = 0
    for [s,e,l, id] in intervals:
        forbidden = {}
        for j in range(i):
            [sp, ep, _, _] = intervals[j]
            if s <= ep:
                forbidden[rank[j]] = True
        r = 0
        while r in forbidden:
            r += 1
        rank.append(r)
        plt.plot([s,e], [r, r], lw=1, color="blue")
        plt.plot([s,s], [0, r], lw=0.3, color="gray")
        plt.plot([e,e], [0, r], lw=0.3, color="gray")
        plt.text((s+e)/2, r+0.05, str(id))
        i += 1
    ax.set_yticks([])
    ax.get_yaxis().set_visible(False)
    plt.savefig("r4.png")
    # plt.show()