require "spec_helper"

RSpec.describe PasvLib do
  it "has a version number" do
    expect(PasvLib::VERSION).not_to be nil
  end

  describe PasvLib::Utils do
    let(:klass) { Class.new.extend PasvLib::Utils }

    describe "#spans_start" do
      it "is true when the query has a non gap before the start" do
        gapped_start = 10
        query = "A--------------------CTGCTG"

        expect(klass.spans_start query, gapped_start).to be true
      end

      it "is true if the query has a non gap at the ROI start" do
        gapped_start = 3
        query = "--ACTG"

        expect(klass.spans_start query, gapped_start).to be true
      end

      it "is false if the query has only gaps before the start" do
        gapped_start = 10
        query = "---------------------CTGCTG"

        expect(klass.spans_start query, gapped_start).to be false
      end
    end

    describe "#spans_end" do
      it "is true when the query has a non gap after the end" do
        gapped_end = 10
        query = "A--------------------CTGCTG"

        expect(klass.spans_end query, gapped_end).to be true
      end

      it "is true if the query has a non gap at the ROI end" do
        gapped_end = 3
        query = "ACT--"

        expect(klass.spans_end query, gapped_end).to be true
      end


      it "is false if the query has only gaps after the end" do
        gapped_end = 10
        query = "ACTG--------------------------"

        expect(klass.spans_end query, gapped_end).to be false
      end
    end

    describe "#pos_to_gapped_pos" do
      it "takes original positions and returns the gapped positions" do
        gapped_key_seq = "A-C-T-G"
        map = { 1 => 1, 2 => 3, 3 => 5, 4 => 7 }

        expect(klass.pos_to_gapped_pos gapped_key_seq).to eq map
      end
    end

    describe "#get_oligo" do
      it "gets the oligo from the query seq given key positions" do
        key_posns = [2, 3]
        gapped_query_seq = "A-C-T-G"
        map = { 1 => 1, 2 => 3, 3 => 5, 4 => 7 }
        oligo = "CT"

        expect(klass.get_oligo gapped_query_seq, key_posns, map).
          to eq oligo
      end
    end

    describe "#get_type" do
      it "builds the type from the oligo and spanning info" do
        oligo = "CT"
        spans = "Yes"

        expect(klass.get_type oligo, spans).to eq "CT_Yes"
      end

      it "leaves off spanning if spanning is NA" do
        oligo = "CT"
        spans = "NA"

        expect(klass.get_type oligo, spans).to eq "CT"
      end
    end
  end
end
