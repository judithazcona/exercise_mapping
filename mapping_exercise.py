import random
from Bio import SeqIO


def random_text(n, alphabet="ACGT"):
    text = ""
    for _ in range(n):
        text += random.choice(alphabet)
    return text + "$"

# --------------------------- TASK 1 ---------------------------
print("\n--- TASK 1 ---\n")

def brute_force(text, p):
    lista = []
    i = 0

    while i + len(p) <= len(text):
        if text[i:i + len(p)] == p:
            lista.append(i)
        i += 1

    return lista

print(brute_force("AAAA", "AA"))

def suffix_array(text):
    sufix = []
    for i in range(len(text)):
        sufix.append((text[i:], i))

    sufix = sorted(sufix)

    lista = []
    for sufijo, i in sufix:
        lista.append(i)

    return lista


# self-checks 
for alphabet in ["ACGT", "AC"]:
    for _ in range(300):
        text = random_text(random.randint(1, 40), alphabet)
        sa = suffix_array(text)
        assert sorted(sa) == list(range(len(text)))
        for i in range(len(sa) - 1):
            assert text[sa[i]:] < text[sa[i + 1]:]


def find_block(text, sa, p):
    n = len(p)

    lo = 0
    hi = len(sa)

    while lo < hi:  # a = first row where the piece is >= p
        mid = (lo + hi) // 2
        start = sa[mid]
        piece = text[start:start + n]
        if piece < p:
            lo = mid + 1
        else:
            hi = mid
    a = lo

    lo = 0
    hi = len(sa)

    while lo < hi:  # b = first row where the piece is > p
        mid = (lo + hi) // 2
        start = sa[mid]
        piece = text[start:start + n]
        if piece <= p:
            lo = mid + 1
        else:
            hi = mid
    b = lo

    return (a, b)


# self-checks 
for alphabet in ["ACGT", "AC"]:
    for _ in range(300):
        text = random_text(random.randint(1, 40), alphabet)
        sa = suffix_array(text)
        p = random_text(random.randint(1, 5), alphabet)[:-1]
        a, b = find_block(text, sa, p)
        assert sorted(sa[a:b]) == brute_force(text, p)

# --------------------------- TASK 2 ---------------------------

print("\n--- TASK 2 ---\n")

def bwt(text, sa):
    result = ""
    for row in range(len(sa)):
        result += text[sa[row] - 1]
    return result

print(bwt("CAT$", [3, 1, 0, 2]))

#self checks
for alphabet in ["ACGT", "AC"]:
    for _ in range(300):
        text = random_text(random.randint(1, 40), alphabet)
        sa = suffix_array(text)
        b = bwt(text, sa)
        assert sorted(b) == sorted(text)
        assert b.count("$") == 1


def first_occurrence(letters):
    counts = {}
    for k in letters:
        if k in counts:
            counts[k] += 1
        else:
            counts[k] = 1

    result = {}
    for X in counts:
        F = 0
        for k in counts:
            if k < X:
                F += counts[k]
        result[X] = F

    return result


print(first_occurrence("TC$A"))

# self checks
for alphabet in ["ACGT", "AC"]:
    for _ in range(300):
        text = random_text(random.randint(1, 40), alphabet)
        fo = first_occurrence(bwt(text, suffix_array(text)))
        assert fo["$"] == 0
        letters = sorted(fo)
        for i in range(len(letters) - 1):
            assert fo[letters[i]] < fo[letters[i + 1]]


def inverse_bwt(bwt):
    seen = {}
    ks = []
    for ch in bwt: 
        seen[ch] = seen.get(ch, 0) + 1
        ks.append(seen[ch])

    fo = first_occurrence(bwt)
    row = 0
    collected = []
    for _ in range(len(bwt) - 1): 
        X = bwt[row]
        k = ks[row]
        collected.append(X)
        row = fo[X] + k - 1
    collected.reverse()
    return "".join(collected) + "$"


print(inverse_bwt("TC$A"))

# self checks
for alphabet in ["ACGT", "AC"]:
    for _ in range(300):
        text = random_text(random.randint(1, 40), alphabet)
        assert inverse_bwt(bwt(text, suffix_array(text))) == text

# think about it
def count_runs(s):
    runs = 0
    for i in range(len(s)):
        if i == 0 or s[i] != s[i - 1]:
            runs += 1
    return runs


ref_record = SeqIO.read("reference.fasta", "fasta")
ref = str(ref_record.seq)

t1 = ref + "$"
t2 = ref + ref + "$"
print("runs text:", count_runs(t1))
print("runs bwt:", count_runs(bwt(t1, suffix_array(t1))))
print("runs doubled text:", count_runs(t2))
print("runs bwt doubled:", count_runs(bwt(t2, suffix_array(t2))))

# --------------------------- TASK 3 ---------------------------

print("\n--- TASK 3 ---\n")

class FMIndex:
    def __init__(self, text, s):
        sa = suffix_array(text)
        n = len(text)
        self.s = s
        self.bwt = bwt(text, sa)
        self.first_occurrence = first_occurrence(self.bwt)

        self.counts = {X: [0] * (n + 1) for X in set(self.bwt)}
        for i, ch in enumerate(self.bwt):
            for X in self.counts:
                self.counts[X][i + 1] = self.counts[X][i]
            self.counts[ch][i + 1] += 1

       
        self.sampled = {}
        for row in range(n):
            if sa[row] % s == 0:
                self.sampled[row] = sa[row]

    def lf(self, row):
        X = self.bwt[row]
        return self.first_occurrence[X] + self.counts[X][row]

    def count(self, p):
        top = 0
        bottom = len(self.bwt)
        for X in reversed(p): 
            if X not in self.first_occurrence: 
                return (0, 0)
            top = self.first_occurrence[X] + self.counts[X][top]
            bottom = self.first_occurrence[X] + self.counts[X][bottom]
            if top == bottom:
                return (top, top)
        return (top, bottom)

    def locate(self, p, _steps=None):
        top, bottom = self.count(p)
        positions = []
        for row in range(top, bottom):
            r = row
            steps = 0
            while r not in self.sampled: 
                r = self.lf(r)
                steps += 1
            if _steps is not None:
                _steps.append(steps)
            positions.append(self.sampled[r] + steps)
        return positions


# examples
fm = FMIndex("CAT$", 2)
print(fm.bwt)
print(fm.counts["T"])
print(fm.sampled)
print(fm.lf(1))
print(fm.count("AT"))
print(fm.locate("AT"))

# self checks
for alphabet in ["ACGT", "AC"]:
    for s in [1, 4, 32]:
        for _ in range(300):
            text = random_text(random.randint(1, 40), alphabet)
            fm = FMIndex(text, s)
            sa = suffix_array(text)
            p = random_text(random.randint(1, 5), alphabet)[:-1]

            top, bottom = fm.count(p)
            steps = []
            loc = fm.locate(p, steps)
            expected = brute_force(text, p)

            assert bottom - top == len(loc) == len(expected)  # 1
            assert sorted(loc) == expected  # 2
            a, b = find_block(text, sa, p)  # 3
            if expected:
                assert (top, bottom) == (a, b)
            else:
                assert top == bottom and a == b
            assert sorted(fm.sampled.values()) == list(range(0, len(text), s))  # 4
            assert all(st < s for st in steps)  # 5
print("Task 3 self-checks OK")

# think about it
fm_ref = FMIndex(ref + "$", 32)
total_ints = sum(len(v) for v in fm_ref.counts.values())
print("lists in counts:", len(fm_ref.counts))
print("integers in counts:", total_ints, "| letters in bwt:", len(fm_ref.bwt))

# --------------------------- TASK 4 ---------------------------

print("\n--- TASK 4 ---\n")


def reverse_complement(seq):
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    result = ""
    for ch in reversed(seq):
        result += comp[ch]
    return result


print(reverse_complement("AACG"))
for _ in range(100):
    x = random_text(random.randint(1, 40))[:-1]
    assert reverse_complement(reverse_complement(x)) == x


def hamming(a, b):
    assert len(a) == len(b)
    d = 0
    for x, y in zip(a, b):
        if x != y:
            d += 1
    return d


print(hamming("ACGT", "AGGA"))


def seeds(read, n):
    L = len(read)
    result = []
    for i in range(n):
        start = i * L // n
        end = (i + 1) * L // n
        result.append((start, read[start:end]))
    return result


print(seeds("ACGTACGTAC", 3))
for n in range(1, 8):  
    r = random_text(random.randint(10, 30))[:-1]
    sd = seeds(r, n)
    assert "".join(piece for _, piece in sd) == r


def map_read(read, fm, ref, d=2, n_seeds=None):
    if n_seeds is None:
        n_seeds = d + 1
    hits = []
    total_candidates = 0
    for strand, seq in (("+", read), ("-", reverse_complement(read))):
        candidates = set() 
        for offset, seed in seeds(seq, n_seeds):
            for p in fm.locate(seed):
                start = p - offset
                if 0 <= start <= len(ref) - len(seq): 
                    candidates.add(start)
        total_candidates += len(candidates)
        for start in candidates: 
            m = hamming(seq, ref[start:start + len(seq)])
            if m <= d:
                hits.append((strand, start, m))
    return sorted(hits), total_candidates


def brute_map(read, ref, d=2):
    hits = []
    for strand, seq in (("+", read), ("-", reverse_complement(read))):
        for start in range(len(ref) - len(seq) + 1):  # every start
            m = hamming(seq, ref[start:start + len(seq)])
            if m <= d:
                hits.append((strand, start, m))
    return sorted(hits)


# map all reads
fm = FMIndex(ref + "$", 32)
reads = [(rec.id, str(rec.seq)) for rec in SeqIO.parse("reads.fasta", "fasta")]
print("reads:", len(reads))

results = {}
total_cand = 0
for name, seq in reads:
    hits, cand = map_read(seq, fm, ref)
    results[name] = hits
    total_cand += cand

mapped = sum(1 for h in results.values() if h)
unmapped = len(reads) - mapped
multiple = sum(1 for h in results.values() if len(h) > 1)
print("mapped:", mapped, "| unmapped:", unmapped, "| more than one hit:", multiple)

for name, seq in reads[:5]:
    for strand, pos, mm in results[name]:
        print(name, strand, pos, mm)
    if not results[name]:
        print(name, "unmapped")

# self checks
L = 100
for name, seq in reads:  # 1
    for strand, pos, mm in results[name]:
        assert 0 <= pos <= len(ref) - L
        s_ = seq if strand == "+" else reverse_complement(seq)
        assert hamming(s_, ref[pos:pos + L]) == mm <= 2

for name, seq in random.sample(reads, 30):  # 2
    assert results[name] == brute_map(seq, ref)

for name, seq in reads:  # 3
    if not results[name]:
        assert brute_map(seq, ref) == []

assert ("+", 0, 0) in map_read(ref[:100], fm, ref)[0]  # 4
assert ("-", len(ref) - 100, 0) in map_read(reverse_complement(ref[-100:]), fm, ref)[0]
print("Task 4 self-checks OK")


baseline = results
print("n_seeds=3 (default): candidates per read =", round(total_cand / len(reads), 2))
for n_seeds in [10, 12, 16]:
    tot = 0
    same = True
    for name, seq in reads:
        hits, cand = map_read(seq, fm, ref, n_seeds=n_seeds)
        tot += cand
        if hits != baseline[name]:
            same = False
    print(f"n_seeds={n_seeds}: seed length ~{100 // n_seeds}, "
          f"candidates per read={tot / len(reads):.2f}, same hits={same}")


kmers = {}
for i in range(len(ref) - 12 + 1):
    kmers[ref[i:i + 12]] = kmers.get(ref[i:i + 12], 0) + 1
print("12-mers that appear more than once:", {k: v for k, v in kmers.items() if v > 1})
print("reads with 2 or more hits:", multiple)

