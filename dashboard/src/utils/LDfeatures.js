import { withLDProvider, withLDConsumer } from "launchdarkly-react-client-sdk";

const user = {
  firstName: 'Shijie',
  lastName: 'Yao',
  key: 'New feature',
  custom: {
    groups: 'beta_testers'
  }
};