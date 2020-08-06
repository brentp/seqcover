var assert = require('assert');
let data = require("./data.js")
let transcript = require("../transcript.js")

let tr = new transcript.Transcript(data.tr_data)

describe('Transcript', function () {
  describe('cdsstart', function () {
    it('should equal cdsstart', function () {
        assert.equal(tr.cdsstart, data.tr_data.cdsstart)
    });
  });
});
