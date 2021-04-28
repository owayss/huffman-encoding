//! # Huffman encoding
//!
//! `huffman-encoding` is a crate to encode information with minimum-redundancy codes
//! using the [Huffman coding](https://en.wikipedia.org/wiki/Huffman_coding)
//! algorithm.
//!
//! ## References
//!
//! * _Hamming, R.R., 1997. Art of doing science and engineering: Learning to learn. CRC Press._
//! * _Huffman, D.A., 1952. A method for the construction of minimum-redundancy codes. Proceedings of the IRE, 40(9), pp.1098-1101._

use bitvec::prelude::*;
use unicode_segmentation::UnicodeSegmentation;

#[derive(Debug, PartialEq)]
pub enum ErrorKind<'a> {
    MissingEncoding(&'a BitSlice),
}

/// An encoded symbol is represented as
/// a [`bitvec::vec::BitVec`](https://docs.rs/bitvec/0.22.3/bitvec/vec/struct.BitVec.html), a contiguous
/// array of bits.
pub type Encoded = BitVec;
/// The encoding's decision tree is represented as a hash map from string keys
/// ([a Unicode grapheme cluster](http://www.unicode.org/reports/tr29/#Grapheme_Cluster_Boundaries))
/// to their [`Encoded`] representation).
pub type Tree = std::collections::HashMap<String, Encoded>;

const ZERO: bool = false;
const ONE: bool = true;

/// Creates and returns an ordered list of pairs of the characters found in the
/// input with their count, ordered in decreasing frequency.
///
/// # Examples
///
/// Basic usage:
///
/// ```
/// let freqs = huffman_encoding::frequencies("huffman");
/// let mut iter = freqs.iter();
///
/// assert_eq!(iter.next(), Some(&("f", 2)));
/// assert_eq!(iter.next(), Some(&("a", 1)));
/// assert_eq!(iter.next(), Some(&("h", 1)));
/// assert_eq!(iter.next(), Some(&("m", 1)));
/// assert_eq!(iter.next(), Some(&("n", 1)));
/// assert_eq!(iter.next(), Some(&("u", 1)));
/// assert_eq!(iter.next(), None);
/// ```
pub fn frequencies<'a>(s: &'a str) -> Vec<(&'a str, usize)> {
    let mut freq = ::std::collections::HashMap::new();

    for g in UnicodeSegmentation::graphemes(s, true) {
        *(freq.entry(g).or_insert(0)) += 1;
    }
    let mut symbols = freq.into_iter().collect::<Vec<(&str, usize)>>();
    symbols.sort_by(|a, b| {
        // The ordering by comparing the characters here in case the
        // frequencies are equal is only to ensure that the algorithm
        // is deterministic, regardless of the order in which the
        // characters appear in the input stream.
        b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0))
    });
    symbols
}

// The forward step in Huffman's algorithm combines the two least
// frequent symbols into one symbol of frequency equal to the sum of
// the frequencies of its parts and inserts it in its new position to maintain
// the list sorted in order of frequencies.
// The process goes on until all symbols are merged into one symbol that
// has the structural hierarchy of the frequencies list.
fn merge<'a>(frequencies: &'a Vec<(&str, usize)>) -> Option<Symbol<'a>> {
    match frequencies.len() {
        0 => None,
        1 => Some(Symbol::Node(
            Some(Box::new(Symbol::Leaf((frequencies[0].0, frequencies[0].1)))),
            None,
        )),
        mut n => {
            let mut heap = ::std::collections::BinaryHeap::with_capacity(n);
            for s in frequencies {
                heap.push(Symbol::Leaf((s.0, s.1)));
            }
            while n > 2 {
                let s = Symbol::Node(
                    Some(Box::new(heap.pop().unwrap())),
                    Some(Box::new(heap.pop().unwrap())),
                );
                heap.push(s);
                n -= 1;
            }
            Some(Symbol::Node(
                Some(Box::new(heap.pop().unwrap())),
                Some(Box::new(heap.pop().unwrap())),
            ))
        }
    }
}

// The backward step in Huffman's algorithm works is given a root symbol with
// all the assembled frequencies and proceeds to unfold hierarchy of symbols:
// at each step, a node is split into its two constituent parts (sub-trees),
// where one is assigned a prefix of `0` and the other a prefix of `1`.
// Whenever a leaf node is encountered, its encoding (the combined prefixes along the path from the root to that node) is emitted.
fn split(
    node: &Symbol,
    prefix: &Encoded,
    encodings: &mut std::collections::HashMap<String, Encoded>,
) {
    match node {
        Symbol::Leaf((s, _)) => {
            encodings.insert(s.to_string(), prefix.to_bitvec());
        }
        Symbol::Node(l, r) => {
            if let Some(l) = l {
                let mut lprefix = prefix.to_bitvec();
                lprefix.push(ZERO);
                split(l, &lprefix, encodings);
            }
            if let Some(r) = r {
                let mut rprefix = prefix.to_bitvec();
                rprefix.push(ONE);
                split(r, &rprefix, encodings);
            }
        }
    }
}

// Constructs the Huffman encoding tree for an input stream s.
fn tree(root: &Option<Symbol>) -> Tree {
    let mut encodings = Tree::new();
    if let Some(root) = root {
        split(&root, &BitVec::new(), &mut encodings);
    }
    encodings
}

/// Encodes an input string using [Huffman's coding](https://en.wikipedia.org/wiki/Huffman_coding)
/// for minimum-redundancy codes of variable length
/// of the encoded input as a bit vecotr and its encoding tree.
///
/// # Examples
///
/// Basic usage:
///
/// ```
/// use huffman_encoding::*;
///
/// let (tree, encoded) = huffman_encoding::encode("baba");
/// {
///     use bitvec::prelude::*;
///     let mut expected_tree = Tree::new();
///     expected_tree.insert("a".to_string(), bitvec![0]);
///     expected_tree.insert("b".to_string(), bitvec![1]);
///     let expected_encoded = bitvec![1, 0, 1, 0];
///     assert_eq!(tree, expected_tree);
///     assert_eq!(encoded, expected_encoded);
/// }
/// ```
pub fn encode(s: &str) -> (Tree, Encoded) {
    let encoding = tree(&merge(&frequencies(s)));
    let mut encoded = Encoded::new();
    for c in UnicodeSegmentation::graphemes(s, true) {
        encoded.extend(encoding.get(c).unwrap());
    }
    (encoding, encoded)
}

/// Decodes an `encoded` string given its `encoding` tree.
///
/// Returns an error if no encoding in the tree was found for a chunk of bits.
///
/// # Examples
///
/// Basic usage:
///
/// ```
/// use huffman_encoding::*;
///
/// let (tree, encoded) = encode("huffman");
/// let decoded = decode(&tree, &encoded).unwrap();
/// assert_eq!(decoded, "huffman");
///
/// ```
///
/// # Errors
///
/// Returns an [`ErrorKind::MissingEncoding`] with the remaining slice of bits
/// that could not be decoded.
///
/// ```
/// use huffman_encoding::*;

/// let (mut tree, encoded) = encode("huffman");
/// let n_encoding = tree.remove("n").unwrap();
/// assert_eq!(
///     decode(&tree, &encoded)
///         .expect_err("should find no encoding for the final symbol 'n')"),
///     ErrorKind::MissingEncoding(&n_encoding)
/// );
/// ```
pub fn decode<'a>(encoding: &Tree, encoded: &'a Encoded) -> Result<String, ErrorKind<'a>> {
    let dict = encoding
        .iter()
        .map(|(sym, enc)| (enc, sym.as_str()))
        .collect::<std::collections::HashMap<&Encoded, &str>>();

    let mut it = encoded.iter().enumerate();
    let mut sym = Encoded::new();
    let mut decoded = Vec::<&str>::new();
    let mut start = 0;
    while let Some((i, c)) = it.next() {
        sym.push(*c);
        if let Some(&s) = dict.get(&sym) {
            decoded.push(s);
            sym.truncate(0);
            start = i + 1;
        }
    }
    if !sym.is_empty() {
        return Err(ErrorKind::MissingEncoding(&encoded[start..]));
    }
    Ok(decoded.join(""))
}

#[derive(Debug, Clone)]
enum Symbol<'a> {
    Leaf((&'a str, usize)),
    Node(Option<Box<Symbol<'a>>>, Option<Box<Symbol<'a>>>),
}

impl Symbol<'_> {
    pub fn frequency(&self) -> usize {
        match self {
            Symbol::Leaf((_, f)) => *f,
            Symbol::Node(None, None) => 0,
            Symbol::Node(Some(l), None) => Symbol::frequency(l),
            Symbol::Node(None, Some(r)) => Symbol::frequency(r),
            Symbol::Node(Some(l), Some(r)) => Symbol::frequency(l) + Symbol::frequency(r),
        }
    }
}
impl Ord for Symbol<'_> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // We flip the order of the arguments as we want the ordering to be
        // decreasing.
        other.frequency().cmp(&self.frequency())
    }
}

impl PartialOrd for Symbol<'_> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for Symbol<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.frequency() == other.frequency()
    }
}
impl Eq for Symbol<'_> {}

#[cfg(test)]
mod tests {
    use crate::*;
    #[test]
    fn frequencies() {
        assert_eq!(crate::frequencies(""), vec![]);
        assert_eq!(crate::frequencies("a"), vec![("a", 1)]);
        assert_eq!(
            crate::frequencies("aaabc"),
            vec![("a", 3), ("b", 1), ("c", 1)]
        );
        assert_eq!(
            crate::frequencies("baaac"),
            vec![("a", 3), ("b", 1), ("c", 1)]
        );
        assert_eq!(
            crate::frequencies("caaab"),
            vec![("a", 3), ("b", 1), ("c", 1)]
        );
        assert_eq!(crate::frequencies("ضَ"), vec![("ضَ", 1)]);
    }

    #[test]
    fn symbol() {
        // Higher frequency symbls come first
        assert!(Symbol::Leaf(("a", 3)) == Symbol::Leaf(("b", 3)));
        assert!(Symbol::Leaf(("a", 3)) < Symbol::Leaf(("b", 2)));
        assert!(Symbol::Leaf(("a", 3)) < Symbol::Leaf(("a", 2)));
        assert!(Symbol::Leaf(("a", 3)) > Symbol::Leaf(("c", 4)));
        assert!(Symbol::Leaf(("a", 3)) > Symbol::Leaf(("a", 4)));
    }

    #[test]
    fn frequnecy() {
        assert_eq!(Symbol::Leaf(("s", 3)).frequency(), 3);
        assert_eq!(
            Symbol::Node(Some(Box::new(Symbol::Leaf(("s", 3)))), None,).frequency(),
            3
        );
        assert_eq!(
            Symbol::Node(None, Some(Box::new(Symbol::Leaf(("s", 3))))).frequency(),
            3
        );
        assert_eq!(
            Symbol::Node(
                Some(Box::new(Symbol::Leaf(("a", 2)))),
                Some(Box::new(Symbol::Node(
                    Some(Box::new(Symbol::Leaf(("b", 1)))),
                    Some(Box::new(Symbol::Leaf(("c", 1))))
                )))
            )
            .frequency(),
            4
        );
    }

    #[test]
    fn merge() {
        assert_eq!(
            crate::merge(&vec![("a", 3), ("b", 1), ("c", 1),]),
            Some(Symbol::Node(
                Some(Box::new(Symbol::Leaf(("a", 3)))),
                Some(Box::new(Symbol::Node(
                    Some(Box::new(Symbol::Node(
                        Some(Box::new(Symbol::Leaf(("c", 1)))),
                        None,
                    ))),
                    Some(Box::new(Symbol::Leaf(("b", 1))))
                )))
            ))
        );
    }

    #[test]
    fn split() {
        let tree = Symbol::Node(
            Some(Box::new(Symbol::Leaf(("a", 3)))),
            Some(Box::new(Symbol::Node(
                Some(Box::new(Symbol::Node(
                    Some(Box::new(Symbol::Leaf(("c", 1)))),
                    None,
                ))),
                Some(Box::new(Symbol::Leaf(("b", 1)))),
            ))),
        );
        {
            let mut expected = std::collections::HashMap::new();
            expected.insert("a".to_string(), bitvec![0]);
            expected.insert("b".to_string(), bitvec![1, 1]);
            expected.insert("c".to_string(), bitvec![1, 0, 0]);
            let mut encodings = std::collections::HashMap::new();
            crate::split(&tree, &bitvec![], &mut encodings);
            assert_eq!(encodings, expected);
        }
    }

    #[test]
    fn tree() {
        {
            assert_eq!(crate::tree(&None), std::collections::HashMap::new());
        }
        {
            let mut encodings = std::collections::HashMap::new();
            encodings.insert("a".to_string(), bitvec![0]);
            let f = crate::frequencies("a");
            let root = crate::merge(&f);
            assert_eq!(crate::tree(&root), encodings);
        }
        {
            let mut encodings = std::collections::HashMap::new();
            encodings.insert("a".to_string(), bitvec![0]);
            encodings.insert("b".to_string(), bitvec![1, 0]);
            encodings.insert("c".to_string(), bitvec![1, 1]);
            let f = crate::frequencies("aaaabbcc");
            let root = crate::merge(&f);
            assert_eq!(crate::tree(&root), encodings);
        }
    }

    #[test]
    fn encode() {
        assert_eq!(crate::encode("").1, bits![]);
        assert_eq!(crate::encode("a").1, bits![0]);
        assert_eq!(crate::encode("ba").1, bits![1, 0]);
        {
            let mut tree = Tree::new();
            tree.insert("a".to_string(), bitvec![0]);
            tree.insert("b".to_string(), bitvec![1]);
            let encoded = bitvec![1, 0, 1, 0];
            assert_eq!(crate::encode("baba"), (tree, encoded));
        }
    }

    #[test]
    fn decode() {
        {
            let mut encoding = Tree::new();
            encoding.insert("h".to_string(), bitvec![0, 1]);
            encoding.insert("e".to_string(), bitvec![0, 0]);
            encoding.insert("l".to_string(), bitvec![1, 1]);
            encoding.insert("o".to_string(), bitvec![1, 0]);
            let encoded = bitvec![0, 1, 0, 0, 1, 1, 1, 1, 1, 0];
            assert_eq!(crate::decode(&encoding, &encoded).unwrap(), "hello");
        }
        {
            let mut encoding = Tree::new();
            encoding.insert("h".to_string(), bitvec![0, 1]);
            encoding.insert("e".to_string(), bitvec![0, 0]);
            encoding.insert("l".to_string(), bitvec![1, 1]);
            encoding.insert("o".to_string(), bitvec![1, 0]);
            let mut encoded = bitvec![0, 1, 0, 0, 1, 1, 1, 1, 1, 0];
            let malformed = bitvec![0];
            encoded.extend(&malformed);
            assert_eq!(
                crate::decode(&encoding, &encoded)
                    .expect_err("should return extra malformed bits [0]"),
                ErrorKind::MissingEncoding(&malformed)
            );
        }
        {
            let mut encoding = Tree::new();
            encoding.insert("a".to_string(), bitvec![0]);
            encoding.insert("b".to_string(), bitvec![1, 0]);
            encoding.insert("c".to_string(), bitvec![1, 1, 0]);
            encoding.insert("d".to_string(), bitvec![1, 1, 1, 0]);
            encoding.insert("e".to_string(), bitvec![1, 1, 1, 1]);
            let mut encoded = bitvec![1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0]; // "edcba"
            let malformed = bitvec![1, 1, 1];
            encoded.extend(&malformed);
            assert_eq!(
                crate::decode(&encoding, &encoded).expect_err(&format!(
                    "should return extra malformed bits {:?}",
                    &malformed
                )),
                ErrorKind::MissingEncoding(&malformed)
            );
        }
    }
}
