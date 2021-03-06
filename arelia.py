#! /usr/bin/env python

import os
import sys
import numpy as np
from numpy.lib.stride_tricks import as_strided

# background probabilites of 20 amino acids from BLOSUM62.
bprob = np.array(
    [
        .078, .051, .041, .052, .024, .034, .059, .083, .025, .062,
        .092, .056, .024, .044, .043, .059, .055, .014, .034, .072
    ]
)

# odd ratio of 20 X 20 amino acid matrix from BLOSUM62.
odd_ratio = np.array([
    [
        3.5339, .5782, .5941, .5424, .8547,
        .7164, .6519, .8959, .5641, .6617,
        .6132, .7555, .6944, .4662, .6559,
        1.3690, .8625, .3663, .4902, .9081
    ],
    [
        .5782, 6.8435, .9565, .6033, .3268,
        1.4418, .8973, .4016, .9412, .3795,
        .5115, 2.1709, .6536, .4011, .4560,
        .7644, .6417, .4202, .5190, .4357
    ],
    [
        .5941, .9565, 8.3879, 1.7355, .4065,
        1.0760, .9095, .8522, 1.3659, .3934,
        .3712, 1.0453, .5081, .4435, .5105,
        1.2815, .9756, .3484, .5022, .4065
    ],
    [
        .5424, .6033, 1.7355, 7.8772, .3205,
        .9050, 1.5971, .5792, .7692, .3722,
        .3135, .8242, .4006, .3497, .5367,
        .9126, .6643, .2747, .3394, .3472
    ],
    [
        .8547, .3268, .4065, .3205, 20.6597,
        .3676, .2825, .4016, .3333, .7392,
        .7246, .3720, .6944, .4735, .3876,
        .7062, .6818, .2976, .3676, .8102
    ],
    [
        .7164, 1.4418, 1.0760, .9050, .3676,
        6.3149, 1.7448, .4961, 1.1765, .4269,
        .5115, 1.6282, .8578, .3342, .5472,
        .9472, .7487, .4202, .6055, .4902
    ],
    [
        .6519, .8973, .9095, 1.5971, .2825,
        1.7448, 4.6251, .3880, .9492, .3280,
        .3685, 1.2409, .4944, .3467, .5518,
        .8618, .6163, .3632, .4487, .4002
    ],
    [
        .8959, .4016, .8522, .5792, .4016,
        .4961, .3880, 5.4870, .4819, .2721,
        .2750, .5379, .3514, .3286, .3923,
        .7760, .4819, .3442, .2835, .3012
    ],
    [
        .5641, .9412, 1.3659, .7692, .3333,
        1.1765, .9492, .4819, 14.8800, .3871,
        .4348, .8571, .6667, .7273, .4651,
        .7458, .5091, .5714, 1.7647, .3333
    ],
    [
        .6617, .3795, .3934, .3722, .7392,
        .4269, .3280, .2721, .3871, 4.7867,
        1.9986, .4608, 1.6801, 1.0997, .3751,
        .4647, .7918, .4608, .6641, 2.6882
    ],
    [
        .6132, .5115, .3712, .3135, .7246,
        .5115, .3685, .2750, .4348, 1.9986,
        4.3833, .4852, 2.2192, 1.3340, .3539,
        .4422, .6522, .5435, .7033, 1.4342
    ],
    [
        .7555, 2.1709, 1.0453, .8242, .3720,
        1.6282, 1.2409, .5379, .8571, .4608,
        .4852, 5.1339, .6696, .3653, .6645,
        .9383, .7468, .3827, .5252, .4712
    ],
    [
        .6944, .6536, .5081, .4006, .6944,
        .8578, .4944, .3514, .6667, 1.6801,
        2.2192, .6696, 6.9444, 1.1364, .3876,
        .6356, .7576, .5952, .7353, 1.3310
    ],
    [
        .4662, .4011, .4435, .3497, .4735,
        .3342, .3467, .3286, .7273, 1.0997,
        1.3340, .3653, 1.1364, 9.4525, .2643,
        .4622, .4959, 1.2987, 2.8075, .8207
    ],
    [
        .6559, .4560, .5105, .5367, .3876,
        .5472, .5518, .3923, .4651, .3751,
        .3539, .6645, .3876, .2643, 10.3299,
        .6701, .5920, .1661, .3420, .3876
    ],
    [
        1.3690, .7644, 1.2815, .9126, .7062,
        .9472, .8618, .7760, .7458, .4647,
        .4422, .9383, .6356, .4622, .6701,
        3.6196, 1.4484, .3632, .4985, .5650
    ],
    [
        .8625, .6417, .9756, .6643, .6818,
        .7487, .6163, .4819, .5091, .7918,
        .6522, .7468, .7576, .4959, .5920,
        1.4484, 4.1322, .3896, .4813, .9091
    ],
    [
        .3663, .4202, .3484, .2747, .2976,
        .4202, .3632, .3442, .5714, .4608,
        .5435, .3827, .5952, 1.2987, .1661,
        .3632, .3896, 33.1633, 1.8908, .3968
    ],
    [
        .4902, .5190, .5022, .3394, .3676,
        .6055, .4487, .2835, 1.7647, .6641,
        .7033, .5252, .7353, 2.8075, .3420,
        .4985, .4813, 1.8908, 8.8235, .6127
    ],
    [
        .9081, .4357, .4065, .3472, .8102,
        .4902, .4002, .3012, .3333, 2.6882,
        1.4342, .4712, 1.3310, .8207, .3876,
        .5650, .9091, .3968, .6127, 3.7809
    ]
])


class Path:
    @staticmethod
    def mkdir_p(path):
        '''
        Make a diretory if not exists.
        '''
        if path:
            import errno
            try:
                os.makedirs(path)
            except OSError as exc:
                if exc.errno == errno.EEXIST and os.path.isdir(path):
                    pass
                else:
                    raise

    @staticmethod
    def listfiles_r(inpath, exts=[], depth_count=1):
        '''
        Iterate over files in a given directory recursively.
        `exts` limits the result with extensions of the files.
        '''
        for filename in os.listdir(inpath):
            filepath = os.path.join(inpath, filename)
            if os.path.isdir(filepath):
                for f in Path.listfiles_r(
                    filepath, exts=exts, depth_count=depth_count+1
                ):
                    yield f
            elif os.path.isfile(filepath):
                if exts and not os.path.splitext(filepath)[1] in exts:
                    continue
                subpath = filename
                stem = os.path.split(filepath)[0]
                for c in range(depth_count-1):
                    stem, base = os.path.split(stem)
                    subpath = os.path.join(base, subpath)
                yield filepath, subpath


class _Seq:
    def __init__(self, tag, seq):
        self.tag = tag
        self.seq = seq


class Seq(list):
    def parse_fasta(self, string):
        tag, seq = '', ''
        for l in string.splitlines():
            if not l:
                continue
            if l[0] == '>':
                if seq:
                    self.add(tag.strip(), seq.strip().replace(' ', ''))
                    seq = ''
                tag = l[1:]
            else:
                seq += l
        self.add(tag.strip(), seq.strip().replace(' ', ''))

    def parse_phylip(self, string):
        blocks = string.split('\n\n')
        i = 0
        for block in blocks:
            if not block:
                continue
            j = 0
            for l in block.splitlines():
                if not l:
                    continue
                if l[0] != ' ':
                    cs = l.split()
                    self.add(cs[0], ''.join(cs[1:]))
                elif i > 0:
                    self[j].seq += l.strip().replace(' ', '')
                j += 1
            i += 1

    def parse_msf(self, string):
        _, alg = string.split('//')
        blocks = alg.strip().split('\n\n')
        tags = []
        tag2seq = {}
        for block in blocks:
            if not block:
                continue
            for l in block.splitlines():
                if not l:
                    continue
                cs = l.split()
                if not cs[0] in tags:
                    tags.append(cs[0])
                    tag2seq[cs[0]] = ''
                tag2seq[cs[0]] += ''.join(cs[1:])
        for tag in tags:
            self.add(tag, tag2seq[tag])

    def make(self, string, infmt='fasta'):
        infmt = infmt.lower()
        if infmt in ('fasta', 'fas', 'mfa', 'fna', 'fsa', 'fa'):
            self.parse_fasta(string)
        elif infmt == 'msf':
            self.parse_msf(string)
        elif infmt in ('phylip', 'phy'):
            self.parse_phylip(string)

    def add(self, tag, seq):
        self.append(_Seq(tag, seq))

    @staticmethod
    def chunks(S, n):
        return [S[c:c+n] for c in range(0, len(S), n)]

    def write_phylip(
        self, oh,
        max_tag_len=10, line_len=60, block_len=10
    ):
        seq_len = len(self[0].seq)
        seq_num = len(self)
        seq_offset = max_tag_len+3
        oh.write(' %s %s\n' % (seq_num, seq_len))
        for i, s in enumerate(range(0, seq_len, line_len)):
            for j in range(seq_num):
                if i:
                    header = ' '*seq_offset
                else:
                    tag = self[j].tag[:max_tag_len]
                    header = tag + ' '*(seq_offset-len(tag))
                seq = ' '.join(
                    self.chunks(self[j].seq[s:s+line_len], block_len)
                )
                oh.write(header+seq+'\n')
            oh.write('\n')

    def write_fasta(self, oh, line_len=60):
        for ob in self:
            oh.write(
                ">%s\n%s\n" % (
                    ob.tag,
                    "\n".join(self.chunks(ob.seq, line_len))
                )
            )

    def write(
        self, oh, outfmt='fasta',
        max_tag_len=10, line_len=60, block_len=10
    ):
        outfmt = outfmt.lower()
        if outfmt in ('phylip','phy'):
            self.write_phylip(
                oh,
                max_tag_len=max_tag_len,
                line_len=line_len,
                block_len=block_len
            )

        elif outfmt == 'tsv':
            for ob in self:
                oh.write(ob.tag+'\t'+ob.seq+'\n')

        elif outfmt == 'csv':
            for ob in self:
                oh.write(ob.tag+','+ob.seq+'\n')

        elif outfmt in ('fasta', 'fas', 'mfa', 'fna'):
            self.write_fasta(oh, line_len=line_len)

    @classmethod
    def parse_string(cls, string, infmt='fasta'):
        o = cls()
        o.make(string, infmt=infmt)
        return o

    @classmethod
    def parse_seqarr(cls, seqarr, tags):
        o = cls()
        nseqlen = seqarr.shape[1]
        if nseqlen:
            for i, seq in enumerate(seqarr.view('S'+str(nseqlen))[:, 0]):
                # decode the frag for python3 compatibility.
                o.add(tags[i], str(seq.decode()))
            return o


# copied from numpy 1.9 to support lower versions
def unique(ar, return_index=False, return_inverse=False, return_counts=False):
    """
    Find the unique elements of an array.
    Returns the sorted unique elements of an array. There are three optional
    outputs in addition to the unique elements: the indices of the input array
    that give the unique values, the indices of the unique array that
    reconstruct the input array, and the number of times each unique value
    comes up in the input array.
    Parameters
    ----------
    ar : array_like
        Input array. This will be flattened if it is not already 1-D.
    return_index : bool, optional
        If True, also return the indices of `ar` that result in the unique
        array.
    return_inverse : bool, optional
        If True, also return the indices of the unique array that can be used
        to reconstruct `ar`.
    return_counts : bool, optional
        If True, also return the number of times each unique value comes up
        in `ar`.
        .. versionadded:: 1.9.0
    Returns
    -------
    unique : ndarray
        The sorted unique values.
    unique_indices : ndarray, optional
        The indices of the first occurrences of the unique values in the
        (flattened) original array. Only provided if `return_index` is True.
    unique_inverse : ndarray, optional
        The indices to reconstruct the (flattened) original array from the
        unique array. Only provided if `return_inverse` is True.
    unique_counts : ndarray, optional
        The number of times each of the unique values comes up in the
        original array. Only provided if `return_counts` is True.
        .. versionadded:: 1.9.0
    """

    ar = np.asanyarray(ar).flatten()

    optional_indices = return_index or return_inverse
    optional_returns = optional_indices or return_counts

    if ar.size == 0:
        if not optional_returns:
            ret = ar
        else:
            ret = (ar,)
            if return_index:
                ret += (np.empty(0, np.bool),)
            if return_inverse:
                ret += (np.empty(0, np.bool),)
            if return_counts:
                ret += (np.empty(0, np.intp),)
        return ret

    if optional_indices:
        perm = ar.argsort(kind='mergesort' if return_index else 'quicksort')
        aux = ar[perm]
    else:
        ar.sort()
        aux = ar
    flag = np.concatenate(([True], aux[1:] != aux[:-1]))

    if not optional_returns:
        ret = aux[flag]
    else:
        ret = (aux[flag],)
        if return_index:
            ret += (perm[flag],)
        if return_inverse:
            iflag = np.cumsum(flag) - 1
            inv_idx = np.empty(ar.shape, dtype=np.intp)
            inv_idx[perm] = iflag
            ret += (inv_idx,)
        if return_counts:
            idx = np.concatenate(np.nonzero(flag) + ([ar.size],))
            ret += (np.diff(idx),)
    return ret


class Profile:
    @staticmethod
    def maparr(arr, K, V):
        '''
        Maps a given numpy array using given keys and values.
        Returns a mapped numpy array.

        Parameters
            arr - a numpy array
            K   - an 1D numpy array (keys)
            V   - an 1D numpy array (values)

        Returns
            the mapped ndarray
        '''
        arg = K.argsort()
        return V[arg][np.searchsorted(K[arg], arr.ravel())].reshape(arr.shape)

    @staticmethod
    def intarr2fracs(intarr, max_int, weights=None):
        '''
        Calculates weighted or unweighted fractions of integers
        in individual columns (axis 0).

        Parameters
            intarr   - a 2D integer numpy array
            max_int  - integers < max_int will be counted
            weights  - none or a 1D numpy array with weights of
                       the same length to intarr axis 0

        Returns
            A 2D numpy array with fractions of integers binned by axis 0
        '''
        seqnum, seqlen = intarr.shape

        frac = np.empty((seqnum, max_int, seqlen), dtype=np.float64)

        if np.any(weights) is not None:
            sw_means = (weights.sum()-weights)

            if np.any(sw_means == 0):
                return Profile.intarr2fracs(intarr, max_int, weights=None)

            for i in range(max_int):
                warr = (intarr == i) * weights[:, None]
                frac[:, i] = (np.sum(warr, axis=0) - warr) / sw_means[:, None]
        else:
            rest_seqnum = seqnum - 1.0
            for i in range(max_int):
                warr = (intarr == i)
                frac[:, i] = (np.sum(warr, axis=0) - warr) / rest_seqnum

        return frac

    @staticmethod
    def intarr2scr(
        intarr, gaparr, max_int,
        scr_types=('s1'), gap_penalty=-5.0, weights=None
    ):
        '''
        Transform a given integer array to profiled scores

        Parameters
            intarr    - a 2D integer numpy array
            gaparr    - a 2D boolean numpy array
            max_int   - integers < max_int will be counted
            scr_types - a tuple with score types('s1', 's2')
            weights   - none or a 1D numpy array with weights of
                       the same length to intarr axis 0

        Returns
            A 2D numpy array with profiled scores
        '''
        seqnum, seqlen = intarr.shape

        fracs = Profile.intarr2fracs(intarr, max_int, weights=weights)

        if 's1' in scr_types:
            p1 = np.empty(intarr.shape)
            p1.fill(gap_penalty)
        if 's2' in scr_types:
            p2 = np.empty(intarr.shape)
            p2.fill(gap_penalty)

        for i in range(fracs.shape[0]):
            IF = ~(gaparr[i])
            nfrac = np.compress(IF, fracs[i], axis=1)
            Ng = 1.0 - nfrac.sum(axis=0)  # weighted gap fraction
            nfrac += Ng * bprob[:, None]
            s = np.log(np.sum(
                np.take(odd_ratio, intarr[i, IF], axis=1) * nfrac, axis=0)
            )
            if 's1' in scr_types:
                p1[i][IF] = s
            if 's2' in scr_types:
                p2[i][IF] = s + gap_penalty * Ng

        r = []
        for each in scr_types:
            if each == 's1':
                r.append((each, p1))
            if each == 's2':
                r.append((each, p2))
        return r

    @staticmethod
    def intarr2weight(intarr):
        '''
        Calculate Position-Specific Sequence Weights from `intarr`.

        Parameters
            intarr - an 2D numpy integer array (< 100)

        Returns
            An 1D numpy float array containing weights for sequences
        '''
        shp = intarr.shape
        us, inv, counts = unique(
            intarr+np.arange(0, shp[1]*100, 100),
            return_inverse=True, return_counts=True
        )
        colind, inv2ind, col_counts = unique(
            (us/100).astype(np.int64),
            return_inverse=True, return_counts=True
        )
        col_w = 1.0 / col_counts
        char_w = 1.0 / counts
        return (col_w[inv2ind] * char_w)[inv].reshape(shp).sum(axis=1)


def sliding_window1d(arr, ws):
    '''
    Return a sliding window over an 1-dimensional numpy array

    Parameters
        arr - an 1-dimensional numpy array
        ws  - an int for window size

    Returns
        an 2D numpy array containing each window from `arr`
    '''

    strides = [arr.strides[0]] * 2
    return as_strided(arr, shape=(arr.shape[0]-ws+1, ws), strides=strides)


def write_csv(fh, X, fmt, delimiter):
    '''
    Write 1D or 2D numpy array

    Parameters
        fh         - writable
        X          - 1- or 2-dimensional numpy array
        fmt        - format
        delimiter  - delimiter

    '''
    shape_len = len(X.shape)
    if shape_len == 2:
        FMT = delimiter.join([fmt]*X.shape[1])
        for x in X:
            fh.write((FMT+'\n') % tuple(x))
    elif shape_len == 1:
        FMT = delimiter.join([fmt]*X.shape[0])
        fh.write((FMT+'\n') % tuple(X))


class ARELIA(dict):
    '''
    Residue and gap characters in the form of 'byte'.
    First 20 characters for amino acids,
    the others for gap and ambiguous symbols.
    '''

    # 20 amino acid residues
    res = 'ARNDCQEGHILKMFPSTWYV'
    gap_eles = 'BZX-*'
    # the number of amino acids
    res_len = len(res)

    # dtype='S' is necessary because python3 string format is unicode.
    ele_bytes = np.array([res+gap_eles], dtype='S').view(np.byte)

    # integer representations of ele_bytes.
    # Ambiguous symbols are treated as a gap.
    ele_ints = np.arange(len(ele_bytes), dtype=np.uint8)
    ele_ints[res_len:] = res_len

    def __init__(self, seq, infmt='fasta', tags=None):
        if infmt == 'bytearr':
            if tags is None or len(tags) != len(seq):
                self.exit('ERROR: The number of tags and seqs must be same.\n')
            self.load_from_bytearr(seq, tags)
        else:
            self.load_from_seqtxt(seq, infmt=infmt)

    def load_from_seqtxt(self, seqtxt, infmt='fasta'):
        tags, seqs = [], []
        seqlen = None
        for e in Seq.parse_string(seqtxt, infmt=infmt):
            cur_seqlen = len(e.seq)
            if cur_seqlen == 0:
                self.exit('ERROR: MSA format (%s) might be wrong.\n' % (infmt))
            if seqlen is not None and seqlen != cur_seqlen:
                self.exit('ERROR: Inconsistent Sequence Length.\n')
            tags.append(e.tag)
            seqs.append([e.seq.upper()])
            seqlen = cur_seqlen

        if len(seqs) < 2:
            self.exit('ERROR: Insufficient sequence number.\n')

        # dtype='S' is necessary because python3 string format is unicode.
        # transform sequences to a 2D integer array with bytes.
        self.load_from_bytearr(np.array(seqs, dtype='S').view(np.byte), tags)

    def load_from_bytearr(self, bytearr, tags):
        self.bytearr = bytearr
        self.tags = np.array(tags)
        # map bytes into integers from 0 to 20
        self.intarr = Profile.maparr(
            self.bytearr, self.ele_bytes, self.ele_ints
        )
        self.alnnum, self.alnlen = self.bytearr.shape
        # make a 2D boolean array that marks gaps
        self.gaparr = self.intarr == self.res_len

    def exit(self, mess=None):
        if mess:
            sys.stderr.write(mess)
        sys.exit()

    @classmethod
    def multi_patch_filter(self, scores, W):
        '''
        Multi-patch filter

        parameters
            scores: 1D numpy array with numeric values
            W: 1D numpy array of window sizes

        Returns
            1D numpy array of filtered values
        '''

        seq_len = scores.shape[0]
        max_w = np.max(W)
        adj = max_w - 1

        nS = np.empty((seq_len+(2*adj)), dtype=np.float64)
        nS.fill(-100.0)  # place very low values at blank positions.
        nS[adj:-adj] = scores
        wS = sliding_window1d(nS, max_w)

        RS = sliding_window1d(wS.sum(axis=-1), max_w).max(axis=-1)/max_w
        for w in W:
            d = max_w - w
            if not d:
                continue
            RS += sliding_window1d(wS[:-d, d:].sum(axis=-1), w).max(axis=-1)/w
        return RS / len(W)

    def get_scr_types(self):
        keys = []
        for key in self:
            if key[0] == 's':
                keys.append(key)
        return keys

    def get_res_scr_types(self):
        keys = []
        for key in self:
            if key[0] == 'r':
                keys.append(key)
        return keys

    def cal_profile_scr(
        self, weighting=True,
        gap_cut_accept=.7, gap_cut_cons=.3, gap_penalty=-5.0,
        scr_types=['s1']
    ):

        gapfrac = self.gaparr.sum(axis=0) / float(self.alnnum)
        self.cons_cols = gapfrac < gap_cut_cons
        self.accept_cols = gapfrac < gap_cut_accept

        weights = None
        if weighting and np.count_nonzero(self.cons_cols) > 0:
            weights = Profile.intarr2weight(
                np.compress(self.cons_cols, self.intarr, axis=1)
                )

        SCR = Profile.intarr2scr(
                    np.compress(self.accept_cols, self.intarr, axis=1),
                    np.compress(self.accept_cols, self.gaparr, axis=1),
                    self.res_len,
                    scr_types=scr_types,
                    gap_penalty=gap_penalty,
                    weights=weights
                )

        for scr_type, score in SCR:
            self[scr_type] = np.empty((self.alnnum, self.alnlen))
            self[scr_type].fill(gap_penalty)
            self[scr_type][:, self.accept_cols] = score

    def cal_res_scr(self, W=[5,10,15,30], gap_penalty=-5.0):
        arg = ~self.gaparr + self.cons_cols  # set patch-fitered positions
        for scr_type in self.get_scr_types():
            rscr_type = 'r'+scr_type
            self[rscr_type] = np.empty(self.bytearr.shape)
            self[rscr_type].fill(gap_penalty)
            for i in range(self.alnnum):
                self[rscr_type][i][arg[i]] = self.multi_patch_filter(
                    self[scr_type][i, arg[i]], W)

    def get_col_scr(self, scr_type):
        if np.any(self.get(scr_type)) is not None:
            RS = self[scr_type].copy()
            RS[self.gaparr] = float('inf')
            return RS.min(axis=0)

    def write_res_scr(self, fh, scr_type, quiet=True):
        if fh:
            write_csv(fh, self[scr_type], '%.2f', ',')
            if not quiet:
                print('saved at '+fh.name)
            fh.close()

    def write_col_scr(self, fh, scr_type, quiet=True):
        if fh:
            CS = self.get_col_scr(scr_type)
            if np.any(CS) is not None:
                write_csv(fh, CS, '%.2f', '\n')
                if not quiet:
                    print('saved at '+fh.name)
            fh.close()

    def mask_res(
        self, RS, cutoff, replacement='-',
        keep_length=False, trim_taxa=False
    ):
        lowarr = RS < cutoff
        nbytearr = self.bytearr.copy()
        nbytearr[np.logical_and(lowarr, ~self.gaparr)] = ord(replacement)

        ngaparr = np.logical_or(self.gaparr, lowarr)
        if keep_length:
            if trim_taxa:
                argt = ~np.all(ngaparr, axis=1)
                return (np.compress(argt, nbytearr, axis=0),
                    np.compress(argt, self.tags))
            else:
                return nbytearr, self.tags
        else:
            argc = ~np.all(ngaparr, axis=0)
            nbytearr = np.compress(argc, nbytearr, axis=1)
            if trim_taxa:
                argt = ~np.all(ngaparr, axis=1)
                return (np.compress(argt, nbytearr, axis=0),
                    np.compress(argt, self.tags))
            else:
                return nbytearr, self.tags

    def mask_col(
        self, CS, cutoff, replacement='-', keep_length=False, trim_taxa=False
    ):
        arg_pass = CS >= cutoff
        if keep_length:
            nbytearr = self.bytearr.copy()
            nbytearr[:, ~arg_pass] = ord(replacement)
            ngaparr = ~arg_pass + self.gaparr
            if trim_taxa:
                argt = ~np.all(ngaparr, axis=1)
                return (np.compress(argt, nbytearr, axis=0),
                    np.compress(argt, self.tags))
            else:
                return nbytearr, self.tags
        else:
            if trim_taxa:
                argt = ~np.all(
                    np.compress(arg_pass, self.gaparr, axis=1),
                    axis=1)
                return (np.compress(
                    argt, np.compress(arg_pass, self.bytearr, axis=1),
                    axis=0),
                    np.compress(arg_pass, self.tags))
            else:
                return np.compress(arg_pass, self.bytearr, axis=1), self.tags

    def get_res_masked_msa(
        self, scr_type, cutoff, replacement='-',
        keep_length=True, trim_taxa=False, quiet=True
    ):
        RS = self.get(scr_type)
        if np.any(RS) is not None:
            return Seq.parse_seqarr(
                *self.mask_res(
                    RS, cutoff, replacement=replacement,
                    keep_length=keep_length, trim_taxa=trim_taxa
                ))

    def get_col_masked_msa(
        self, scr_type, cutoff, replacement='-',
        keep_length=False, trim_taxa=False, quiet=True
    ):
        CS = self.get_col_scr(scr_type)
        if np.any(CS) is not None:
            return Seq.parse_seqarr(
                *self.mask_col(
                    CS, cutoff, replacement=replacement,
                    keep_length=keep_length, trim_taxa=trim_taxa
                    ))

    def write_res_masked_msa(
        self, fh, scr_type, cutoff,
        replacement='-', outfmt='fasta',
        keep_length=True, trim_taxa=False, quiet=True
    ):
        if fh:
            msa = self.get_res_masked_msa(
                scr_type, cutoff, replacement=replacement,
                keep_length=keep_length, trim_taxa=trim_taxa
                )
            if msa:
                msa.write(fh, outfmt=outfmt)
                if not quiet:
                    print('saved at '+fh.name)
                self._res_msa_saved = True
            else:
                self._res_msa_saved = False
            fh.close()

    def write_col_masked_msa(
        self, fh, scr_type, cutoff,
        replacement='-', outfmt='fasta',
        keep_length=False, trim_taxa=False, quiet=True
    ):
        if fh:
            msa = self.get_col_masked_msa(
                scr_type, cutoff,
                replacement=replacement, keep_length=keep_length,
                trim_taxa=trim_taxa, quiet=True
                )
            if msa:
                msa.write(fh, outfmt=outfmt)
                if not quiet:
                    print('saved at '+fh.name)
                self._col_msa_saved = True
            else:
                self._col_msa_saved = False
            fh.close()

    @classmethod
    def go(
            cls, seqtxt, W=[5, 10, 15, 30], weighting=False,
            cutoff=0.3, gap_cut_accept=.7, gap_cut_cons=.3,
            gap_penalty=-5.0, scr_type_res='s1', scr_type_col='s2',
            msa_col=None, msa_res=None, scr_col=None, scr_res=None,
            infmt='fasta', outfmt='fasta', replacement='-',
            keep_length=False, trim_taxa=False,
            write=False, quiet=False,
     ):
        if (
            msa_col is None and msa_res is None and
            scr_col is None and scr_res is None
        ):

            # set msa_res as stdout when no print-out options are given.
            msa_res = sys.stdout
            quiet = True

        scr_types = []
        if msa_col or scr_col:
            scr_types.append(scr_type_col)
        if msa_res or scr_res:
            scr_types.append(scr_type_res)

        arelia = cls(seqtxt, infmt=infmt)
        arelia.cal_profile_scr(
            weighting=weighting,
            gap_cut_accept=gap_cut_accept,
            gap_cut_cons=gap_cut_cons,
            gap_penalty=gap_penalty,
            scr_types=scr_types
        )

        if 's1' in arelia or 's2' in arelia:
            arelia.cal_res_scr(W=W, gap_penalty=gap_penalty)
            if write:
                arelia.write_res_masked_msa(
                    msa_res, 'r'+scr_type_res, cutoff,
                    replacement=replacement, outfmt=outfmt,
                    keep_length=keep_length, trim_taxa=trim_taxa,
                    quiet=quiet
                )
                arelia.write_res_scr(
                    scr_res, 'r'+scr_type_res,
                    quiet=quiet
                )
                arelia.write_col_masked_msa(
                    msa_col, 'r'+scr_type_col, cutoff,
                    replacement=replacement, outfmt=outfmt,
                    keep_length=keep_length, trim_taxa=trim_taxa, quiet=quiet
                )
                arelia.write_col_scr(scr_col, 'r'+scr_type_col, quiet=quiet)

        return arelia


if __name__ == '__main__':
    import argparse

    # arguments
    parser = argparse.ArgumentParser(
        description='''examples:
  -get a residue-masked MSA.
    arelia.py MSA_FILE_IN > MSA_FILE_OUT
   or
    arelia.py MSA_FILE_IN -msa_res MSA_FILE_OUT

  -get a residue-masked MSA and residue reliability scores.
    arelia.py MSA_FILE_IN -msa_res MSA_FILE_OUT -scr_res SCORE_FILE_OUT

  -process all SA files in a directory recursively.
    arelia.py MSA_DIR_IN -msa_res MSA_DIR_OUT -scr_res SCORE_DIR_OUT

  -set a gap penalty, window sizes, and cutoff (residues with scores < cutoff will be masked).
    arelia.py MSA_FILE_IN -msa_res MSA_FILE_OUT -gap -5.0 -W 5 10 15 30 -cutoff 0.3

  -set input and output MSA formats.
    arelia.py MSA_DIR_IN -msa_res MSA_DIR_OUT -infmt fasta -outfmt phylip

''',
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        'input', metavar='PATH',
        help='MSA file or dir containing MSA files.'
    )
    parser.add_argument(
        '-msa_res', metavar='PATH',
        help='file or dir for saving residue masked MSA(s).'
    )
    parser.add_argument(
        '-msa_col', metavar='PATH',
        help='file or dir for saving column masked MSA(s).'
    )
    parser.add_argument(
        '-scr_res', metavar='PATH',
        help='file or dir for saving residue scores.'
    )
    parser.add_argument(
        '-scr_col', metavar='PATH',
        help='file or dir for saving column scores.'
    )
    parser.add_argument(
        '-cutoff', metavar='FLOAT',
        help='cutoff for residue|column masking. default=0.3',
        type=float, default=0.3
    )
    parser.add_argument(
        '-gap', metavar='FLOAT', help='gap penalty. default=-5',
        type=float, default=-5.0
    )
    parser.add_argument(
        '-scr_type_col', metavar='s1|s2',
        help='score type for column scores. default=s2',
        default='s1', choices=['s1', 's2']
    )
    parser.add_argument(
        '-W', metavar='INTEGERS',
        help='window size combination. default=5 10 15 30',
        nargs='+', type=int, default=[5, 10, 15, 30]
    )
    parser.add_argument(
        '-infmt', metavar='fasta|phylip|msf',
        help='MSA input format. default=fasta', default='fas',
        choices=['fasta', 'fas', 'mfa', 'phylip', 'msf']
    )
    parser.add_argument(
        '-outfmt', metavar='fasta|phylip|csv|tsv',
        help='MSA output format. default=fasta', default='fas',
        choices=['fasta', 'fas', 'mfa', 'phylip', 'csv', 'tsv']
    )
    parser.add_argument(
        '-replacement', metavar='S1',
        help='by which unreliable residues will be replaced. default=-',
        default='-'
    )
    parser.add_argument(
        '--keep_length',
        help='keep the column length of an input MSA.',
        action='store_true'
    )
    parser.add_argument(
        '--trim_taxa',
        help='remove taxa with no remaining residues.',
        action='store_true'
    )
    parser.add_argument('--quiet', help='be quiet.', action='store_true')
    parser.add_argument(
        '--version', help='show version.',
        action='version', version='0.3.2'
    )
    args = parser.parse_args()

    def open_p(filename, mode):
        if filename:
            Path.mkdir_p(os.path.split(filename)[0])
            return open(filename, mode)

    def join_p(*args):
        try:
            return os.path.join(*args)
        except:
            return None

    if os.path.isfile(args.input):

        ARELIA.go(
            open(args.input, 'r').read(), args.W, weighting=True,
            cutoff=args.cutoff, gap_cut_accept=.7,
            gap_cut_cons=.3, gap_penalty=args.gap,
            scr_type_col=args.scr_type_col,
            msa_col=open_p(args.msa_col, 'w'),
            msa_res=open_p(args.msa_res, 'w'),
            scr_col=open_p(args.scr_col, 'w'),
            scr_res=open_p(args.scr_res, 'w'),
            infmt=args.infmt,
            outfmt=args.outfmt,
            keep_length=args.keep_length,
            trim_taxa=args.trim_taxa,
            replacement=args.replacement,
            quiet=args.quiet,
            write=True
        )

    elif os.path.isdir(args.input):
        if not args.scr_res and not args.scr_col \
             and not args.msa_res and not args.msa_col:
            parser.error(
                '''%s is a dir. at least one output dir is required.
                   \'arelia.py -h\' will show examples''' % args.input)

        for filename, subfilename in Path.listfiles_r(args.input):
            stem, fname = os.path.split(subfilename)

            ARELIA.go(
                open(filename).read(), args.W, weighting=True,
                cutoff=args.cutoff, gap_cut_accept=.7,
                gap_cut_cons=.3, gap_penalty=args.gap,
                scr_type_col=args.scr_type_col,
                scr_res=open_p(
                    join_p(args.scr_res, stem, fname+'.res.scr'), 'w'
                ),
                scr_col=open_p(
                    join_p(args.scr_col, stem, fname+'.col.scr'), 'w'
                ),
                msa_res=open_p(
                    join_p(args.msa_res, stem, fname+'.res.'+args.outfmt), 'w'
                ),
                msa_col=open_p(
                    join_p(args.msa_col, stem, fname+'.col.'+args.outfmt), 'w'
                ),
                infmt=args.infmt,
                outfmt=args.outfmt,
                keep_length=args.keep_length,
                trim_taxa=args.trim_taxa,
                replacement=args.replacement,
                quiet=args.quiet,
                write=True
            )
