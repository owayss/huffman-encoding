use proptest::prelude::*;
use unicode_segmentation::*;

proptest! {
    #[test]
    fn frequencies(input in any::<String>()) {
        let freqs = huffman_encoding::frequencies(input.as_str());
        let graphemes = UnicodeSegmentation::graphemes(input.as_str(), true).collect::<Vec<&str>>();
        // The sum of the frequencies of all the characters is equal to the
        // length of the input.
        assert_eq!(freqs.iter().fold(0, |acc, ch| acc + ch.1), graphemes.len());
        let graphemes = graphemes.into_iter().collect::<::std::collections::HashSet::<&str>>();
        // The cadinality of the frequencies vector is equal to that of the set
        // of the characters of the input.
        assert_eq!(freqs.len(), graphemes.len());
        // All the elements of the set of the characters of the input
        // are present in the frequencies vector.
        graphemes.iter().for_each(|&g| assert!(freqs.iter().find(|x| x.0 == g).is_some()));
        // The frequencies list is sorted in descending order.
        (1..freqs.len()).for_each(|i| assert!(freqs[i].1 <= freqs[i-1].1))
    }

    #[test]
    fn encode(input in any::<String>()) {
        let (encoding, _encoded) = huffman_encoding::encode(input.as_str());
        // All the characters in the input have an encoding.
        UnicodeSegmentation::graphemes(input.as_str(), true)
            .collect::<::std::collections::HashSet::<&str>>()
            .iter()
            .for_each(|&g| assert!(!encoding.get(g).unwrap().is_empty()));
        // Kraft's inequality holds:
        // https://en.wikipedia.org/wiki/Kraft%E2%80%93McMillan_inequality
        let krafts_sum: f64 = encoding.values().fold(0.0, |acc, enc| acc + (1.0 / (1 << enc.len()) as f64));
        match encoding.len() {
            0 => {},
            1 => assert!(krafts_sum == 0.5),
            _ => assert!(krafts_sum == 1.0),
        }
        // The codes are instantaneously deocdable if no symbol is a prefix to
        // another.
        encoding.iter()
            .for_each(|(k1, v1)| encoding.iter().for_each(|(k2, v2)| assert!(!v2.starts_with(v1) || k1 == k2)));
    }

    #[test]
    fn e2e(input in any::<String>()) {
        let (encoding, encoded) = huffman_encoding::encode(&input);
        assert_eq!(huffman_encoding::decode(&encoding, &encoded).unwrap(), input);
    }
}
