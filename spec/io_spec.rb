require "spec_helper"

RSpec.describe PasvLib::Io do
  let(:klass) { Class.new { extend PasvLib::Io } }

  describe "#read_queries" do
    context "bad queries" do
      context "queries have gaps" do
        let(:infile) do
          # Just reuse the refs file here.
          File.join SpecHelper::TEST_FILE_DIR, "refs_with_gaps.fa"
        end

        it "raises an error" do
          expect {
            klass.read_queries infile
          }.to raise_error PasvLib::ParseError
        end
      end
    end

    context "good queries" do
      let(:infile) do
        File.join SpecHelper::TEST_FILE_DIR, "small_refs.fa"
      end

      let(:expected_queries) do
        {
          "pasv_query___s1 apple pie" => "actg",
          "pasv_query___s2 is really"  => "aaccttgg",
          "pasv_query___s3 quite good"  => "aaaccctttggg",
        }
      end

      it "returns refs with special headers" do
        actual_queries = klass.read_queries infile

        expect(actual_queries).to eq expected_queries
      end
    end
  end

  describe "#read_refs" do
    context "bad refs" do
      context "refs have gaps" do
        let(:infile) do
          File.join SpecHelper::TEST_FILE_DIR, "refs_with_gaps.fa"
        end

        it "raises an error" do
          expect {
            klass.read_refs infile
          }.to raise_error PasvLib::ParseError
        end
      end
    end

    context "good refs" do
      context "single reference" do
        let(:infile) do
          File.join SpecHelper::TEST_FILE_DIR, "single_ref.fa"
        end

        let(:expected_refs) do
          {
            "first_pasv_ref" => "actg",
          }
        end

        it "returns refs with special headers" do
          actual_refs = klass.read_refs infile

          expect(actual_refs).to eq expected_refs
        end
      end

      context "multiple references" do
        let(:infile) do
          File.join SpecHelper::TEST_FILE_DIR, "small_refs.fa"
        end

        let (:expected_refs) do
          {
            "first_pasv_ref" => "actg",
            "pasv_ref___s2"  => "aaccttgg",
            "pasv_ref___s3"  => "aaaccctttggg",
          }
        end

        it "returns refs with special headers" do
          actual_refs = klass.read_refs infile

          expect(actual_refs).to eq expected_refs
        end
      end
    end
  end
end