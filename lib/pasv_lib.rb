require "pasv_lib/error"
require "pasv_lib/alignment"
require "pasv_lib/version"
require File.join __dir__, "..", "vendor", "systemu"

module PasvLib
  module CoreExtensions
    module Time
      def date_and_time fmt="%F %T.%L"
        Object::Time.now.strftime fmt
      end

      def time_it title="", logger=nil, run: true
        if run
          t = Object::Time.now

          yield

          time = Object::Time.now - t

          if title == ""
            msg = "Finished in #{time} seconds"
          else
            msg = "#{title} finished in #{time} seconds"
          end

          if logger
            logger.info msg
          else
            $stderr.puts msg
          end
        end
      end
    end

    module Process
      include CoreExtensions::Time

      def run_it *a, &b
        exit_status, stdout, stderr = SystemUniversal::Run.systemu *a, &b


        puts stdout unless stdout.empty?
        unless stderr.empty? || stderr.include?("Insecure world writable dir")
          $stderr.puts stderr
        end

        exit_status
      end

      def run_it! *a, &b
        exit_status = self.run_it *a, &b

        # Sometimes, exited? is not true and there will be no exit
        # status. Success should catch all failures.
        AbortIf.abort_unless exit_status.success?,
                             "Command failed with status " +
                             "'#{exit_status.to_s}' " +
                             "when running '#{a.inspect}', " +
                             "'#{b.inspect}'"

        exit_status
      end

      # Examples
      #
      # Process.extend CoreExtensions::Process
      # Time.extend CoreExtensions::Time
      #
      # Process.run_and_time_it! "Saying hello",
      #                          %Q{echo "hello world"}
      #
      # Process.run_and_time_it! "This will raise SystemExit",
      #                          "ls arstoeiarntoairnt" do
      #   puts "i like pie"
      # end
      def run_and_time_it! title="",
                           cmd="",
                           logger=AbortIf::logger,
                           &b

        AbortIf.logger.debug { "Running: #{cmd}" }

        time_it title, logger do
          run_it! cmd, &b
        end
      end
    end
  end

  module Utils
    # @note Takes 1-based coordinates.
    def spans_start query_seq, gapped_start
      !query_seq[0..gapped_start-1].tr("-", "").empty?
    end

    # @note Takes 1-based coordinates.
    def spans_end query_seq, gapped_end
      !query_seq[gapped_end-1..query_seq.length-1].tr("-", "").empty?
    end

    # @note Returns a hash that accepts and returns 1-based coordinates.
    def pos_to_gapped_pos gapped_key_seq
      pos_to_gapped_pos = {}
      nongap_idx = 0

      gapped_key_seq.each_char.with_index do |char, gapped_idx|
        unless char == "-"
          pos_to_gapped_pos[nongap_idx + 1] = (gapped_idx + 1)

          nongap_idx += 1
        end
      end

      pos_to_gapped_pos
    end

    # @note The key posns are 1-based and non-gapped.  The query seq
    #   is gapped.
    def get_oligo gapped_query_seq, key_posns, pos_to_gapped_pos
      gapped_key_posns = key_posns.map do |pos|
        pos_to_gapped_pos[pos]
      end

      gapped_key_posns.map do |pos|
        idx = pos - 1

        gapped_query_seq[idx]
      end.join.upcase
    end

    def get_type oligo, spans
      if spans == "NA"
        oligo
      else
        "#{oligo}_#{spans}"
      end
    end
  end
end
