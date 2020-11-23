/*
 * Copyright (c) 2020, Oracle and/or its affiliates. All rights reserved.
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 *
 * This code is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 2 only, as
 * published by the Free Software Foundation.
 *
 * This code is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * version 2 for more details (a copy is included in the LICENSE file that
 * accompanied this code).
 *
 * You should have received a copy of the GNU General Public License version
 * 2 along with this work; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * Please contact Oracle, 500 Oracle Parkway, Redwood Shores, CA 94065 USA
 * or visit www.oracle.com if you need additional information or have any
 * questions.
 *
 */

#ifndef SHARE_JFR_RECORDER_SERVICE_JFREVENTTHROTTLER_HPP
#define SHARE_JFR_RECORDER_SERVICE_JFREVENTTHROTTLER_HPP

#include "jfrfiles/jfrEventIds.hpp"
#include "jfr/support/jfrAdaptiveSampler.hpp"

class Semaphore;

class JfrEventThrottler : public JfrAdaptiveSampler {
  friend class JfrRecorder;
  friend class JfrRecorderService;
  friend class JfrEventThrottlerSettingsLocker;
 private:
  JfrSamplerParams _last_params;
  Semaphore* _semaphore;
  double _probability;
  intptr_t _rate_per_second;
  intptr_t _nth_selection;
  JfrEventId _event_id;
  bool _disabled;
  bool _update;

  const JfrSamplerParams& update_params(const JfrSamplerWindow* expired);

  static bool create();
  static void destroy();
  static void clear();
  const JfrSamplerParams& next_window_params(const JfrSamplerWindow* expired) override;
 public:
  JfrEventThrottler(JfrEventId event_id);
  virtual ~JfrEventThrottler();
  bool initialize() override;
  bool configure(intptr_t rate_per_second, double probability, intptr_t nth_selection);
  static JfrEventThrottler* for_event(JfrEventId event_id);
  static bool accept(JfrEventId event_id, int64_t timestamp);
};

#endif // SHARE_JFR_RECORDER_SERVICE_JFREVENTTHROTTLER_HPP