fn main() {}

mod BLAKE2B {
    // Word size in bits, 64 bits for BLAKE2b.
    // https://datatracker.ietf.org/doc/html/rfc7693#section-2.1
    const W: u32 = 64;

    // Number of rounds in compress function F
    // https://datatracker.ietf.org/doc/html/rfc7693#section-3.2
    const ROUND_COUNT: usize = 12;

    // Initializationn vector from SHA-512.
    // https://datatracker.ietf.org/doc/html/rfc7693#section-2.6
    // https://datatracker.ietf.org/doc/html/rfc6234#section-6.3
    const IV: [u64; 8] = [
        0x6A09E667F3BCC908,
        0xBB67AE8584CAA73B,
        0x3C6EF372FE94F82B,
        0xA54FF53A5F1D36F1,
        0x510E527FADE682D1,
        0x9B05688C2B3E6C1F,
        0x1F83D9ABFB41BD6B,
        0x5BE0CD19137E2179,
    ];

    // Message word schedule permutations SIGMA for each round.
    // https://datatracker.ietf.org/doc/html/rfc7693#section-2.7
    const SIGMA: [[usize; 16]; 12] = [
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
        [14, 10, 4, 8, 9, 15, 13, 6, 1, 12, 0, 2, 11, 7, 5, 3],
        [11, 8, 12, 0, 5, 2, 15, 13, 10, 14, 3, 6, 7, 1, 9, 4],
        [7, 9, 3, 1, 13, 12, 11, 14, 2, 6, 5, 10, 4, 0, 15, 8],
        [9, 0, 5, 7, 2, 4, 10, 15, 14, 1, 11, 12, 6, 8, 3, 13],
        [2, 12, 6, 10, 0, 11, 8, 3, 4, 13, 7, 5, 15, 14, 1, 9],
        [12, 5, 1, 15, 14, 13, 4, 10, 0, 7, 6, 3, 9, 2, 8, 11],
        [13, 11, 7, 14, 12, 1, 3, 9, 5, 0, 15, 4, 8, 6, 2, 10],
        [6, 15, 14, 9, 11, 3, 0, 8, 12, 2, 13, 7, 1, 4, 10, 5],
        [10, 2, 8, 4, 7, 6, 1, 5, 15, 11, 9, 14, 3, 12, 13, 0],
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
        [14, 10, 4, 8, 9, 15, 13, 6, 1, 12, 0, 2, 11, 7, 5, 3],
    ];

    // Rotation mixing constants
    // https://datatracker.ietf.org/doc/html/rfc7693#section-2.1
    const R: [u64; 4] = [32, 24, 16, 63];

    // Circular shift of "x" by "n" bits.
    // https://en.wikipedia.org/wiki/Circular_shift
    // https://datatracker.ietf.org/doc/html/rfc7693#section-2.3
    // https://datatracker.ietf.org/doc/html/rfc7693#appendix-C.2
    fn rotate_right_64(x: u64, n: u64) -> u64 {
        (x >> n) ^ (x << (64 - n))
    }

    // Mixing function G, mixes two input words "x" and "y" into four words
    // "a", "b", "c" and "d" in the working vector "v".
    // https://datatracker.ietf.org/doc/html/rfc7693#section-3.1
    fn mix(v: &mut [u64; 16], a: usize, b: usize, c: usize, d: usize, x: u64, y: u64) {
        v[a] = (v[a] + v[b] + x) % 2_u64.pow(W);
        v[d] = rotate_right_64(v[d] ^ v[a], R[0]);
        v[c] = (v[c] + v[d]) % 2_u64.pow(W);
        v[b] = rotate_right_64(v[b] ^ v[c], R[1]);
        v[a] = (v[a] + v[b] + y) % 2_u64.pow(W);
        v[d] = rotate_right_64(v[d] ^ v[a], R[2]);
        v[c] = (v[c] + v[d] + y) % 2_u64.pow(W);
        v[b] = rotate_right_64(v[b] ^ v[c], R[3]);
    }

    // Compression function F takes the state vector "h" and message
    // block vector "m" (last block padded with 0s), 2w-bit offset
    // counter "t" and final block indicator "f".
    // https://datatracker.ietf.org/doc/html/rfc7693#section-3.2
    fn compress(h: &mut [u64; 8], m: &[u64; 16], t: u64, f: bool) {
        // initialize local work vector
        let mut v: [u64; 16] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        v[..7].copy_from_slice(h); // first half from state
        v[8..].copy_from_slice(&IV); // second half from IV

        v[12] = v[12] ^ (t % 2_u64.pow(W));
        v[13] = v[13] ^ (t >> W);

        if f == true {
            v[14] = v[14] ^ 0xFFFFFFFFFFFFFFFF // invert all bits
        }

        // cryptographic mixing
        for i in 0..(ROUND_COUNT - 1) {
            // message word selection permutation for this round
            let s: [usize; 16] = SIGMA[i % 10].clone();

            mix(&mut v, 0, 4, 8, 12, m[s[0]], m[s[1]]);
            mix(&mut v, 1, 5, 9, 13, m[s[2]], m[s[3]]);
            mix(&mut v, 2, 6, 10, 14, m[s[4]], m[s[5]]);
            mix(&mut v, 3, 7, 11, 15, m[s[6]], m[s[7]]);

            mix(&mut v, 0, 5, 10, 15, m[s[8]], m[s[9]]);
            mix(&mut v, 1, 6, 11, 12, m[s[10]], m[s[11]]);
            mix(&mut v, 2, 7, 8, 13, m[s[12]], m[s[13]]);
            mix(&mut v, 3, 4, 9, 14, m[s[14]], m[s[15]]);
        }

        for i in 0..7 {
            h[i] = h[i] ^ v[i] ^ v[i + 8]; // xor the two halves
        }
    }
}
