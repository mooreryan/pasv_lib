require "spec_helper"

require "blosum"

RSpec.describe PasvLib::Alignment do
  let(:klass) { Class.new { extend PasvLib::Alignment }}

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
        seqs = %W[ABCDE]
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
        seqs = ["A B - C", "ab-c"]
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

        expect{klass.alignment_columns seqs}.to raise_error PasvLib::Error
      end
    end
  end

  describe "#similarity_score" do
    let(:max_score) { 1.0 }
    let(:min_score) { 0.0 }

    context "with no sequences" do
      it "raises an error" do
        expect{klass.similarity_score []}.to raise_error PasvLib::Error
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
        seq = "AAAAA"
        seq_gapped = "AA-AA"

        two_seq_score = klass.similarity_score [seq, seq_gapped]
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

        actual_points = 5 + 7 + 2 + 6 -  2 -  2
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
          expect{klass.similarity_score ["AAA", "---"]}.to raise_error PasvLib::Error
        end
      end

      context "alternative chars as gap characters" do
        it "treats '.' as a gap" do
          s1 = "AR.N"
          s2 = "A-DN"
          s3 = ".RNC".downcase

          actual_points = 5 + 7 + 2 + 6 -  2 -  2
          max_points    = 5 + 7 + 7 + 6 + 12 + 12

          similarity_score = klass.similarity_score [s1, s2, s3], Blosum::BLOSUM45

          expect(similarity_score).to eq actual_points / max_points.to_f
        end
      end
    end
  end
end