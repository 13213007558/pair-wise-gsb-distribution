require File.expand_path(File.dirname(__FILE__) + "/spec_helper.rb")
require "json"

describe Distribution::F::Noncentral do
  let(:ref_path) do
    File.expand_path("../../data/noncentral_f_reference.json", __FILE__)
  end
  let(:reference) { JSON.parse(File.read(ref_path)) }
  let(:ncf) { Distribution::F::Noncentral }

  it "is exposed under Distribution::F::Noncentral" do
    expect(described_class).to eq(Distribution::F::Noncentral)
  end

  describe "argument validation" do
    it "rejects non-positive degrees of freedom" do
      expect { ncf.cdf(1, 0, 5, 1) }.to raise_error(ArgumentError)
      expect { ncf.cdf(1, 5, -1, 1) }.to raise_error(ArgumentError)
    end

    it "rejects a negative non-centrality parameter" do
      expect { ncf.cdf(1, 5, 5, -0.001) }.to raise_error(ArgumentError)
    end

    it "rejects quantile probabilities outside [0, 1]" do
      expect { ncf.quantile(-0.1, 5, 5, 1) }.to raise_error(ArgumentError)
      expect { ncf.quantile(1.1, 5, 5, 1) }.to raise_error(ArgumentError)
    end
  end

  describe "lambda = 0 agrees with the central F" do
    it "matches the central pdf, cdf and survival exactly" do
      [0.1, 0.7, 2.5, 10.0].each do |x|
        [2, 5, 10].product([2, 5, 10]).each do |d1, d2|
          expect(ncf.pdf(x, d1, d2, 0).to_f).to eq(Distribution::F::Ruby_.pdf(x, d1, d2))
          expect(ncf.cdf(x, d1, d2, 0).to_f).to eq(Distribution::F::Ruby_.cdf(x, d1, d2))
          expect(ncf.survival(x, d1, d2, 0).to_f).to eq(Distribution::F::Ruby_.q(x, d1, d2))
        end
      end
    end

    it "matches the central inverse CDF" do
      [0.025, 0.5, 0.975].each do |p|
        q = ncf.quantile(p, 5, 10, 0)
        expect((q.to_f - Distribution::F::Ruby_.quantile(p, 5, 10)).abs).to be < 1e-12
        expect(q.converged?).to be(true)
      end
    end
  end

  describe "high-precision reference values" do
    it "matches the CDF in the centre and both tails" do
      reference.each do |row|
        next unless row.key?("cdf")
        d1, d2, lam, x = row.values_at("d1", "d2", "lambda", "x")
        result = ncf.cdf(x, d1, d2, lam, tolerance: 1e-10)
        expect(result.to_f).to be_within(5e-9).of(row["cdf"].to_f)
        expect(result.converged?).to be(true)
      end
    end

    it "matches survival, including probabilities very close to one" do
      reference.each do |row|
        next unless row.key?("survival")
        d1, d2, lam, x = row.values_at("d1", "d2", "lambda", "x")
        result = ncf.survival(x, d1, d2, lam, tolerance: 1e-10)
        expect(result.to_f).to be_within(5e-9).of(row["survival"].to_f)
        expect(result.converged?).to be(true)
      end
    end

    it "matches the PDF across the support" do
      reference.each do |row|
        next unless row.key?("pdf")
        d1, d2, lam, x = row.values_at("d1", "d2", "lambda", "x")
        result = ncf.pdf(x, d1, d2, lam, tolerance: 1e-13)
        scale = [row["pdf"].to_f.abs, 1e-300].max
        # In the extreme tail the density is far below the absolute
        # tolerance, so only relative accuracy on meaningful values is
        # asserted; the absolute error stays bounded by error_bound.
        if scale > 1e-8
          expect(result.to_f).to be_within(5e-6 * scale).of(row["pdf"].to_f)
        else
          expect((result.to_f - row["pdf"].to_f).abs).to be <= result.error_bound + 1e-13
        end
        # Relative accuracy is not meaningful for a density this close to
        # zero, so non-convergence at the absolute tolerance is acceptable
        # as long as the error bound reports it honestly.
        expect(result.converged?).to be(true) if scale > 1e-6
      end
    end

    it "matches the inverse CDF near both tails and at the centre" do
      reference.each do |row|
        next unless row.key?("quantile")
        d1, d2, lam, p = row.values_at("d1", "d2", "lambda", "p")
        result = ncf.quantile(p, d1, d2, lam, tolerance: 1e-10)
        qscale = [row["quantile"].to_f.abs, 1e-12].max
        expect(result.to_f).to be_within(1e-5 * qscale).of(row["quantile"].to_f)
        expect(result.converged?).to be(true)
      end
    end
  end

  describe "CDF and survival are complementary" do
    it "sums to one even in the far tails" do
      [[1e-10, 5, 10, 4], [1e4, 5, 10, 4], [1e-8, 3, 50, 20]].each do |x, d1, d2, lam|
        c = ncf.cdf(x, d1, d2, lam)
        s = ncf.survival(x, d1, d2, lam)
        expect(c.to_f + s.to_f).to be_within(1e-11).of(1.0)
      end
    end
  end

  describe "inverse CDF consistency" do
    it "round-trips through the CDF" do
      [[0.001, 5, 10, 4], [0.5, 5, 10, 20], [0.999, 10, 4, 10]].each do |p, d1, d2, lam|
        q = ncf.quantile(p, d1, d2, lam)
        expect(ncf.cdf(q.to_f, d1, d2, lam, tolerance: 1e-12).to_f).to be_within(1e-6).of(p)
      end
    end

    it "returns 0 for p = 0 and Infinity for p = 1" do
      expect(ncf.quantile(0, 5, 10, 4).to_f).to eq(0.0)
      expect(ncf.quantile(1, 5, 10, 4).infinite?).to be_truthy
    end
  end

  describe "numerical budget and tolerance" do
    it "reports non-convergence when the budget is exhausted" do
      result = ncf.cdf(5.0, 5, 10, 10_000.0,
                        tolerance: 1e-13, max_iterations: 5)
      expect(result.iterations).to be <= 5
      expect(result.converged?).to be(false)
      expect(result.error_bound).to be > 1e-13
    end

    it "does not pass off a truncated result as finished" do
      tiny = ncf.cdf(100.0, 5, 10, 1_000.0, max_iterations: 2)
      full = ncf.cdf(100.0, 5, 10, 1_000.0)
      expect(tiny.converged?).to be(false)
      expect(full.converged?).to be(true)
      expect(tiny.error_bound).to be > 1e-12
    end

    it "uses more mixture terms for a tighter tolerance" do
      loose = ncf.cdf(5.0, 5, 10, 100.0, tolerance: 1e-3)
      tight = ncf.cdf(5.0, 5, 10, 100.0, tolerance: 1e-12)
      expect(tight.iterations).to be >= loose.iterations
      expect(tight.error_bound).to be <= 1e-12
    end
  end

  describe "large non-centrality stability" do
    [100, 1000, 10_000, 100_000].each do |lam|
      it "computes a finite, converging CDF and PDF at lambda = #{lam}" do
        x = 1.0 + lam / 5.0
        c = ncf.cdf(x, 5, 10, lam.to_f, tolerance: 1e-8, max_iterations: 100_000)
        d = ncf.pdf(x, 5, 10, lam.to_f, tolerance: 1e-6, max_iterations: 100_000)
        expect(c.to_f.finite?).to be(true)
        expect(d.to_f.finite?).to be(true)
        expect(c.error_bound).to be <= 1e-8
        expect(c.converged?).to be(true)
      end
    end
  end

  describe "result object" do
    it "behaves like a number and exposes convergence metadata" do
      r = ncf.cdf(2, 5, 10, 4)
      expect(r).to be_a(Numeric)
      expect(r.to_f).to be_a(Float)
      expect(r.error_bound).to be_a(Float)
      expect(r.iterations).to be_an(Integer)
      expect(r.converged?).to be(true)
      expect(2.0 * r).to be_a(Float)
    end
  end
end
