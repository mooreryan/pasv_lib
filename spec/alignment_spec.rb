require "spec_helper"

RSpec.describe PasvLib::Alignment do
  let(:klass) { Class.new { extend PasvLib::Alignment }}

  describe "#get_alignment_columns" do
    it "returns an ary for each column in the alignment" do
      seqs = %w[ABC- AB-D A-CD -BCD]

      expected = [
        %w[A A A -],
        %w[B B - B],
        %w[C - C C],
        %w[- D D D]
      ]

      actual = klass.get_alignment_columns seqs

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

        actual = klass.get_alignment_columns seqs

        expect(actual).to eq expected
      end
    end
  end
end