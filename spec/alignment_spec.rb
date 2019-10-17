require "spec_helper"

require "blosum"

# Calculate the geometric_index by hand up here so it doesn't clutter up the tests.
def expected_geometric_index
  s1 = 'A-N--' # 01011
  s2 = 'AR-D-' # 00101
  s3 = 'AR-D-' # 00101
  s4 = 'A--D-' # 01101

  gene_sequences = [s1, s2, s3, s4]

  # For the geometric homogeneity index.  This looks like a lot, but it is good to double check that the by-hand calculation matches the code, and to list it out explicitly so it is clearer how the algorithm works.
  num_seqs                    = gene_sequences.count
  num_aln_cols                 = gene_sequences.first.length
  max_similarities_per_aln_col = num_seqs
  max_similarities_per_seq     = num_aln_cols
  num_comparisons_per_aln_col  = num_aln_cols - 1
  num_comparisons_per_seq      = num_seqs - 1

  # First do the pairwise comparisons....

  # 1v2, 1v3, 1v4, (1v5 for by residue)
  aln_col1_similarities = [2, 1, 3, 0]
  seq1_similarities = [2, 2, 3]

  # 2v1, 2v3, 2v4, (2v5 for by residue)
  aln_col2_similarities = [2, 1, 3, 2]
  seq2_similarities = [2, 5, 4]

  # 3v1, 3v2, 3v4, (3v5 for by residue)
  aln_col3_similarities = [1, 1, 0, 3]
  seq3_similarities = [2, 5, 4]

  # 4v1, 4v2, 4v3, (4v5 for by residue)
  aln_col4_similarities = [3, 3, 0, 1]
  seq4_similarities = [3, 4, 4]

  # 5v1, 5v2, 5v3, 5v4
  aln_col5_similarities = [0, 2, 3, 1]

  # Then each column has a similarity score w.r.t. the other columns.
  aln_col1_similarity_score = aln_col1_similarities.sum.to_f / max_similarities_per_aln_col / num_comparisons_per_aln_col
  aln_col2_similarity_score = aln_col2_similarities.sum.to_f / max_similarities_per_aln_col / num_comparisons_per_aln_col
  aln_col3_similarity_score = aln_col3_similarities.sum.to_f / max_similarities_per_aln_col / num_comparisons_per_aln_col
  aln_col4_similarity_score = aln_col4_similarities.sum.to_f / max_similarities_per_aln_col / num_comparisons_per_aln_col
  aln_col5_similarity_score = aln_col5_similarities.sum.to_f / max_similarities_per_aln_col / num_comparisons_per_aln_col

  # Also, each seq has a similarity score w.r.t. the other sequences.
  seq1_similarity_score = seq1_similarities.sum.to_f / max_similarities_per_seq / num_comparisons_per_seq
  seq2_similarity_score = seq2_similarities.sum.to_f / max_similarities_per_seq / num_comparisons_per_seq
  seq3_similarity_score = seq3_similarities.sum.to_f / max_similarities_per_seq / num_comparisons_per_seq
  seq4_similarity_score = seq4_similarities.sum.to_f / max_similarities_per_seq / num_comparisons_per_seq

  # The quick geo score is the mean of all alignment column similarity scores.
  by_residue_score = [
    aln_col1_similarity_score,
    aln_col2_similarity_score,
    aln_col3_similarity_score,
    aln_col4_similarity_score,
    aln_col5_similarity_score
  ].sum.to_f / num_aln_cols

  # The by sequence similarity score is the mean of all sequence similarity scores.
  by_sequence_score = [
    seq1_similarity_score,
    seq2_similarity_score,
    seq3_similarity_score,
    seq4_similarity_score
  ].sum.to_f / num_seqs

  # Finally the full geometric similarity score is the average of the quick score (by residue) and the by sequence score.
  (by_residue_score + by_sequence_score) / 2
end

RSpec.describe PasvLib::Alignment do
  let(:klass) { Class.new { extend PasvLib::Alignment } }

  describe "#adjust_scoring_matrix" do
    it "rescales so the min. score is zero" do
      scoring_matrix = {
        "A" => { "A" => 4, "C" => -2 },
        "C" => { "C" => 6, "A" => -3 }
      }

      expected = {
        "A" => { "A" => 7, "C" => 1 },
        "C" => { "C" => 9, "A" => 0 }
      }

      expect(klass.adjust_scoring_matrix scoring_matrix).to eq expected
    end

    context "min value is > 0" do
      it "returns exact copy of matrix if min > 0" do
        scoring_matrix = {
          "A" => { "A" => 4, "C" => 2 },
          "C" => { "C" => 6, "A" => 3 }
        }

        expected = {
          "A" => { "A" => 4, "C" => 2 },
          "C" => { "C" => 6, "A" => 3 }
        }

        expect(klass.adjust_scoring_matrix scoring_matrix).to eq expected
      end

      it "returns a deep copy/new object, not the same thing" do
        scoring_matrix = {
          "A" => { "A" => 4, "C" => 2 },
          "C" => { "C" => 6, "A" => 3 }
        }

        new_matrix = klass.adjust_scoring_matrix scoring_matrix

        scoring_matrix["A"]["A"] = 40000

        expect(new_matrix["A"]["A"]).to eq 4
      end
    end
  end

  describe "#alignment_columns" do
    it "returns an ary for each column in the alignment" do
      seqs = %w[ABC- AB-D A-CD -BCD]

      expected = [
        %w[A A A -],
        %w[B B - B],
        %w[C - C C],
        %w[- D D D]
      ]

      actual = klass.alignment_columns seqs

      expect(actual).to eq expected
    end

    context "input is a single sequence" do
      it "returns a seq_len X 1 array of arrays" do
        seqs     = %W[ABCDE]
        expected = [
          %w[A],
          %w[B],
          %w[C],
          %w[D],
          %w[E]
        ]

        actual = klass.alignment_columns seqs

        expect(actual).to eq expected
      end
    end

    context "alignments with spaces" do
      it "ignores spaces in alignment" do
        seqs     = ["A B - C", "ab-c"]
        expected = [
          %w[A a],
          %w[B b],
          %w[- -],
          %w[C c]
        ]

        actual = klass.alignment_columns seqs

        expect(actual).to eq expected
      end
    end

    context "seqs are not all the same length" do
      it "raises an error" do
        seqs = %w[abcde abc]

        expect { klass.alignment_columns seqs }.to raise_error PasvLib::Error
      end
    end
  end

  describe "#gap?" do
    it "is truthy if residue is a gap" do
      expect(klass.gap? "-").to be_truthy
      expect(klass.gap? ".").to be_truthy
    end

    it "is falsey if residue is a gap" do
      expect(klass.gap? "A").to be_falsey
      expect(klass.gap? "r").to be_falsey
    end
  end

  describe "#geometric_index" do
    it "calculates the geometric index" do
      s1 = 'A-N--' # 01011
      s2 = 'AR-D-' # 00101
      s3 = 'AR-D-' # 00101
      s4 = 'A--D-' # 01101

      expect(klass.geometric_index [s1, s2, s3, s4]).to be_within(1e-5).of expected_geometric_index
    end
  end

  describe "#similarity_score" do
    let(:max_score) { 1.0 }
    let(:min_score) { 0.0 }

    context "with no sequences" do
      it "raises an error" do
        expect { klass.similarity_score [] }.to raise_error PasvLib::Error
      end
    end

    context "with a single sequence" do
      it "has the maximum score" do
        expect(klass.similarity_score ["AAAAA"]).to eq max_score
      end
    end

    context "aligments with identical sequences" do
      it "returns a perfect score" do
        seq = "AAAAA"
        aln = [seq, seq, seq]

        result = klass.similarity_score aln

        expect(result).to eq max_score
      end
    end

    context "handling aligments with gaps" do
      it "ignores gaps while scoring" do
        alns = [
          ["A-A", "A-A", "A-A"],
          ["AAA", "AAA", "A-A"],
          ["-AA", "A-A", "AA-"],
        ]

        alns.each do |aln|
          expect(klass.similarity_score aln).to eq max_score
        end
      end

      it "gives the same score even as you add more gapped seqs" do
        seq        = "AAAAA"
        seq_gapped = "AA-AA"

        two_seq_score   = klass.similarity_score [seq, seq_gapped]
        three_seq_score = klass.similarity_score [seq, seq_gapped, seq_gapped]

        expect(two_seq_score).to eq three_seq_score
      end
    end

    context "similar sequences" do
      it "returns the similarity score" do
        s1 = "AR-N"
        s2 = "A-DN"
        s3 = "-RNC".downcase

        # blosum45
        #
        # score (-2/7) means -2 points out of 7
        #
        # col1 AA (5/5), A- (0/0), A- (0/0)
        # col2 R- (0/0), RR (7/7), -R (0/0)
        # col3 -D (0/0), -N (0/0), DN (2/7) max for D is 7, max for N is 6, so 7 is the max
        # col4 NN (6/6), NC (-2/12), NC (-2/12) max for N is 6, max for C is 12, so 12 is the max

        actual_points = 5 + 7 + 2 + 6 - 2 - 2
        max_points    = 5 + 7 + 7 + 6 + 12 + 12

        similarity_score = klass.similarity_score [s1, s2, s3], Blosum::BLOSUM45

        expect(similarity_score).to eq actual_points / max_points.to_f
      end

      it "can take alternate similarity matrices" do
        scoring_matrix = {
          "A" => { "A" => 4, "C" => 2 },
          "C" => { "C" => 6, "A" => 2 }
        }

        s1 = "AC"
        s2 = "AA"

        expected = (4 + 2) / (4 + 6).to_f

        expect(klass.similarity_score [s1, s2], scoring_matrix).to eq expected
      end
    end

    context "weird edge cases" do
      context "when no comparisons can be made" do
        it "raises an error" do
          expect { klass.similarity_score ["AAA", "---"] }.to raise_error PasvLib::Error
        end
      end

      context "alternative chars as gap characters" do
        it "treats '.' as a gap" do
          s1 = "AR.N"
          s2 = "A-DN"
          s3 = ".RNC".downcase

          actual_points = 5 + 7 + 2 + 6 - 2 - 2
          max_points    = 5 + 7 + 7 + 6 + 12 + 12

          similarity_score = klass.similarity_score [s1, s2, s3], Blosum::BLOSUM45

          expect(similarity_score).to eq actual_points / max_points.to_f
        end
      end

      context "residues that don't exist in the scoring matrix" do
        it "raises a useful error" do
          s1 = "ABCDE"

          scoring_matrix = {
            "A" => { "A" => 4, "C" => 2 },
            "C" => { "C" => 6, "A" => 2 }
          }

          expect { klass.similarity_score [s1, s1], scoring_matrix }.to raise_error PasvLib::Error
        end
      end
    end
  end

  describe "#to_binary_matrix" do
    it "converts an alignment to a binary matrix" do
      s1 = "A-B-"
      s2 = "AA-C"
      s3 = "AB-D"

      aln = [s1, s2, s3]

      expected = [
        [0, 1, 0, 1],
        [0, 0, 1, 0],
        [0, 0, 1, 0],
      ]

      expect(klass.to_binary_matrix aln).to eq expected
    end
  end
end