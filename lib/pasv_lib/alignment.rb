require "blosum"
require "set"

module PasvLib
  module Alignment
    GAP_CHARS = Set.new %w[- .]

    # If the overall min of the scoring matrix is < 0, then this scales it so that the overall min becomes zero.
    #
    # @param scoring_matrix The scoring matrix to rescale.  E.g., Blosum::BLOSUM62.
    #
    # @return A new hash scaled to zero, or a deep copy of the old one if the overall min >= 0.
    def adjust_scoring_matrix scoring_matrix
      overal_min = scoring_matrix.values.map(&:values).flatten.min

      scaling_value = overal_min < 0 ? overal_min.abs : 0

      # We want a deep copy to prevent things from getting weird while using the new hash later.
      new_matrix = {}
      scoring_matrix.each do |residue, scores|
        new_matrix[residue] = {}
        scores.each do |other_residue, score|
          new_matrix[residue][other_residue] = score + scaling_value
        end
      end

      new_matrix
    end

    # Get the columns of an aligned set of sequences.
    #
    # Any spaces in the alignment are ignored (e.g., like those spaces NCBI puts in their alignments sometimes).  Gap characters are '.' or '-'
    #
    # @param seqs [Array<String>] an array of sequences, normally aligned sequences
    #
    # @return [Array<Array<String>>] an array of alignment columns (which are arrays of single residue strings)
    #
    # @raise PasvLib::Error if seqs.empty?
    # @raise PasvLib::Error if no comparisons can be made, e.g., all gaps
    # @raise PasvLib::Error if not all sequences are the same length
    #
    #
    # @example
    #   klass.alignment_columns ["AA-", "A-A"] #=> [["A", "A"], ["A", "-"], ["-", "A"]]
    def alignment_columns seqs
      seqs_no_spaces = seqs.map { |seq| seq.tr " ", "" }
      len = seqs_no_spaces.first.length

      seqs_no_spaces.map do |seq|
        unless seq.length == len
          raise PasvLib::Error, "Aligned seqs must be the same length"
        end

        seq.chars
      end.transpose
    end

    # Calculate the geometric index for an alignment.
    #
    # Basically, you change all residues to 0 and all gaps to 1.  Then you take the permutations of the sequences and then the residues and XOR each of the permutations.  Then you add it up and take averages and you'll get the score.  That is a pretty bad explanation, so see http://merenlab.org/2016/11/08/pangenomics-v2/#geometric-homogeneity-index for more information.
    #
    # @param aln [Array<String>] aligned sequenecs
    # @param by [String] either "sequence" or "residue".  Controls whether to do the calculation by sequence, or by residue.
    #
    # @return [Float] a score between 0 and 1, with 1 being very homogeneous and 1 being very heterogeneous.
    #
    # @note The original Anvi'o code uses a clever bitshifting scheme to avoid storing array of arrays.  It may also speed up the XOR part as you're doing fewer XORs, but I'm not positive about that.
    def geometric_score aln, by
      binary_aln = to_binary_matrix aln

      if by == "residue"
        binary_aln = binary_aln.transpose
      end

      num_rows = binary_aln.length
      max_differences_per_row = binary_aln.first.length
      num_comparisions_per_row = num_rows - 1

      diff_score = binary_aln.permutation(2).map do |(row1, row2)|
        row1.zip(row2).map do |elem1, elem2|
          elem1 ^ elem2
        end.sum / max_differences_per_row.to_f
      end.sum / num_comparisions_per_row.to_f / num_rows

      1 - diff_score
    end

    # A wrapper for #geometric_score that takes the average of the by sequence and by residue scores for an alignment.
    #
    # @param aln [Array<String>] aligned sequenecs
    #
    # @return [Float] a score between 0 and 1, with 1 being very homogeneous and 1 being very heterogeneous.
    def geometric_index aln
      by_seq_score = geometric_score aln, "sequence"
      by_residue_score = geometric_score aln, "residue"

      (by_seq_score + by_residue_score) / 2.0
    end

    # @note Technically you could get a negative score if you don't have enough high scoring residue pairs to keep the total score above zero.  If this is the case, you're alignment probably isn't very good.
    def similarity_score seqs, scoring_matrix = Blosum::BLOSUM62
      raise PasvLib::Error if seqs.empty?
      return 1.0 if seqs.count == 1

      aln_cols = alignment_columns seqs

      actual_points = 0
      max_points    = 0

      aln_cols.each do |residues|
        residues.map(&:upcase).combination(2).each do |r1, r2|
          unless GAP_CHARS.include?(r1) || GAP_CHARS.include?(r2)
            r1_max_score = scoring_matrix[r1].values.max
            r2_max_score = scoring_matrix[r2].values.max
            pair_max     = [r1_max_score, r2_max_score].max

            # TODO check that residues exist in the scoring matrix.
            actual_points += scoring_matrix[r1][r2]
            max_points    += pair_max
          end
        end
      end

      if max_points.zero?
        raise PasvLib::Error, "Something went wrong and max_points was 0.  " \
                              "Maybe one of your sequences is all gaps?"

      end

      actual_points / max_points.to_f
    end

    # Convert an aligment to a binary matrix where 1 represents gaps and 0 represents residues.
    #
    # @param aln [Array<String>] an array of (aligned) sequences
    #
    # @return [Array<Array<integer>>] an array of arrays.  Each row is a sequence.
    def to_binary_matrix aln
      aln.map do |seq|
        seq.chars.map do |char|
          # Gaps turn to 1, residues turn to 0.
          GAP_CHARS.include?(char) ? 1 : 0
        end
      end
    end
  end
end