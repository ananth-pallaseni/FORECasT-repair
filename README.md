# FORECasT-Repair: A repair-aware predictor of Cas9 edeiting outcomes

## Installation
The `forecast_repair` package is dependent on the [self-target](https://github.com/felicityallen/SelfTarget) package, so this must be installed first:

```bash
# Install python dependencies
git clone https://github.com/felicityallen/SelfTarget.git
cd SelfTarget/selftarget_pyutils
pip install -e .
cd ../indel_prediction
pip install -e .

# Compile indel generator
cd indel_analysis/indelmap
cmake . -DINDELMAP_OUTPUT_DIR=/usr/local/bin
make && make install
export INDELGENTARGET_EXE=/usr/local/bin/indelgentarget
```
Make sure that the `INDELGENTARGET_EXE` variable points to the `indelgentarget` binary. If you want the binary to be installed into a different directory, change `/usr/local/bin` in the `cmake` command to point to the directory of choice. 

Once `self-target` has been installed, the `forecast_repair` package is installed as follows:
```bash
git clone https://github.com/ananth-pallaseni/FORECasT-repair.git
cd FORECasT-repair
pip install -e .
```

## Basic usage
```python
import forecast_repair

# First input is a target sequence containing a 20nt spacer, an NGG PAM, and at least 10nt of context on either side of the spacer
target_sequence = 'AATCGCCTAGGTTTCGGCCTCAGTCCAACCACGAGAAGTACAGGGCTTTTACCTTTGCAATCCGGATGGATGTGACGGA'

# Second input is the 0-indexed location of the PAM  inside the target sequence
pam = 42

# Final input is the knockout background in which to make the prediiction
ko = 'Lig4'

# Predictions take the form of a table with the columns: Mutation, Inserted Sequence, and Prediction.
predictions = forecast_repair.predict(target_sequence, pam, ko)

```

