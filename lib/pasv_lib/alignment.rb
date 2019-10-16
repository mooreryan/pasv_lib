module PasvLib
  module Alignment
    # @param seqs [Array<String>] an array of sequences, normally aligned sequences
    #
    # @return [Array<Array<String>>] an array of alignment columns (which are arrays of single residue strings)
    #
    # @example
    #   klass.get_alignment_columns ["AA-", "A-A"] #=> [["A", "A"], ["A", "-"], ["-", "A"]]
    def get_alignment_columns seqs
      seqs.map(&:chars).transpose
    end
  end
end