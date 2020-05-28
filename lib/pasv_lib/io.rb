require "parse_fasta"

module PasvLib
  module Io
    def read_refs fname
      refs = {}

      ParseFasta::SeqFile.open(fname).each_record do |rec|
        if has_gaps? rec.seq
          raise PasvLib::ParseError,
                "Record '#{rec.header}' had gaps!  Did you accidentally " \
                "provide aligned sequences?"
        end

        if refs.count.zero?
          head = "first_pasv_ref"
        else
          head = "pasv_ref___#{rec.id}"
        end

        refs[head] = rec.seq
      end

      refs
    end

    def read_queries fname
      queries = {}

      ParseFasta::SeqFile.open(fname).each_record do |rec|
        if has_gaps? rec.seq
          raise PasvLib::ParseError,
                "Record '#{rec.header}' had gaps!  Did you accidentally " \
                "provide aligned sequences?"
        end

        header          = "pasv_query___#{rec.header}"
        queries[header] = rec.seq
      end

      queries
    end

    private

    def has_gaps? seq
      seq.include?("-") || seq.include?(".")
    end
  end
end